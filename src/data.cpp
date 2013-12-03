// Copyright (c) 2013 Kamil Slowikowski
// See LICENSE for GPLv3 license.

#include "snpspec.h"

// Include functions for controlling threads through OpenMP.
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 0
#define omp_set_num_threads() 0
#endif

// Main function that executes all of the intermediate steps.
snpspec::snpspec(
    std::string user_snpset_file,
    std::string gene_matrix_file,
    std::string gene_intervals_file,
    std::string snp_intervals_file,
    std::string null_snps_file,
    std::string condition_file,
    std::string out_folder,
    ulong slop,
    int threads,
    ulong null_snpset_replicates,
    ulong min_observations,
    ulong max_iterations 
)
{
    std::cout << "# SNPspec v0.1\n\n";

    write_args(
        user_snpset_file,
        gene_matrix_file,
        gene_intervals_file,
        snp_intervals_file,
        null_snps_file,
        condition_file,
        out_folder,
        slop,
        threads,
        null_snpset_replicates,
        min_observations,
        max_iterations,
        std::cout
    );

    // Read names of null SNPs that will be sampled to create random or
    // matched SNP sets.
    std::cout << timestamp() << " # Reading files ..." << std::endl;
    read_names(null_snps_file, _null_snp_names);

    // Optional condition file to condition on specified columns in the
    // gene matrix.
    if (condition_file.length() > 0) {
        read_names(condition_file, _condition_names);
    }

    // Read SNP names and intervals.
    read_bed_intervals(snp_intervals_file, _snp_intervals);

    // Read the gene matrix.
    read_gct(gene_matrix_file, _row_names, _col_names, _gene_matrix);

    // Read the gene intervals but only keep the ones listed in the GCT.
    read_bed_interval_tree(
        gene_intervals_file,
        _row_names,
        _gene_interval_tree
    );

    std::cout << timestamp() << " # done." << std::endl;

    // Report names from the conditions file that are absent from the
    // gene matrix file.
    report_missing_conditions();

    // Drop all SNP intervals except those in the null set.
    //drop_snp_intervals();

    // Check if the matrix is binary by reading the first column.
    if (is_binary(_gene_matrix.col(0))) {
        // Let the user know we detected it.
        std::cout << timestamp() << " # Expression is binary." << std::endl;
        // Cache these values ahead of time.
        _binary_sums = _gene_matrix.colwise().sum();
        _binary_probs = _binary_sums / _gene_matrix.rows();
        _binary_gene_matrix = true;
    } else {
        _binary_gene_matrix = false;

        // Condition the matrix on the specified columns.
        condition(_gene_matrix, _condition_names);

        // Normalize the matrix.
        _gene_matrix =
            _gene_matrix.array().colwise() /
            _gene_matrix.rowwise().norm().eval().array();

        // Reverse percentile rank each column of the matrix.
        // So, a small value like 0.02 means the given gene is highly specific
        // to the column. A large value means the gene is non-specific.
        for (int i = 0; i < _gene_matrix.cols(); i++) {
            _gene_matrix.col(i) =
                rankdata(_gene_matrix.col(i)) / _gene_matrix.rows();
        }
    }

    // 1. Find a geneset for each SNP by querying the gene interval tree.
    // 2. Bin genesets by size. (This will be used to generate SNP sets.)
    const ulong MAX_GENES = 10;
    bin_genesets(slop, MAX_GENES);

    // Check for enrichment of each column in parallel.
    // Ensure that a valid number of threads is used.
    threads = clamp(threads, 1, cpu_count());
    omp_set_num_threads(threads);

    int n_random_snps = 0;

    if (file_exists(user_snpset_file)) {
        read_names(user_snpset_file, _user_snp_names);
    } else {
        random_snps(user_snpset_file, _user_snp_names, slop);
        n_random_snps = _user_snp_names.size();
    }

    mkpath(out_folder);

    std::ofstream args (out_folder + "/args.txt");
    write_args(
        user_snpset_file,
        gene_matrix_file,
        gene_intervals_file,
        snp_intervals_file,
        null_snps_file,
        condition_file,
        out_folder,
        slop,
        threads,
        null_snpset_replicates,
        min_observations,
        max_iterations,
        args
    );
    args.close();

    // Overlap the user's SNP intervals with the gene intervals. Record
    // the SNPs that are not present in the --snp-intervals file. Also
    // record the gene sets and their sizes.
    overlap_genes(
        _user_snp_names,
        _user_absent_snp_names,
        _user_genesets,
        _user_geneset_sizes,
        slop
    );

    // Merge SNPs that share genes or have overlapping genes.
    merge_user_snps(
        _user_snp_names,
        _user_genesets,
        _user_geneset_sizes
    );

    // Report the genes overlapping the user's SNPs.
    report_user_snp_genes(out_folder + "/snp_genes.txt");

    for (auto & size : _user_geneset_sizes) {
        if (size > MAX_GENES) {
            size = MAX_GENES;
        }
    }

    std::cout << timestamp()
              << " # On each iteration, we'll test "
              << _user_geneset_sizes.size()
              << " gene sets from these bins:" << std::endl;
    // Report how many genesets exist of each size.
    for (auto item : _geneset_bins) {
        int n_items = std::count(
            _user_geneset_sizes.begin(),
            _user_geneset_sizes.end(),
            item.first
        );
        if (n_items > 0) {
            std::cout << timestamp()
                      << " # " << setw(3) << n_items
                      << " gene sets with size ";
            if (item.first == MAX_GENES) {
                std::cout << ">= "<< setw(2) << item.first;
            } else {
                std::cout << setw(2) << item.first;
            }
            std::cout << " from a pool of size " << item.second.size()
                      << std::endl;
        }
    }

    std::cout << timestamp()
              << " # Computing up to "
              << scientific << double(max_iterations)
              << " iterations for each column with "
              << fixed << threads << " threads ...\n"
              << std::flush;

    if (null_snpset_replicates > 0) {
        std::cout << timestamp()
                  << " # Computing " 
                  << scientific << null_snpset_replicates
                  << " null SNP sets ...\n"
                  << std::flush;

        for (ulong replicate = 0;
                replicate < null_snpset_replicates; replicate++) {
            // The user specified something like "random20" so let's
            // generate a totally random list of SNPs without matching.
            if (n_random_snps > 0) {
                // Calculate p-values for random null gene sets.
                calculate_pvalues(
                    out_folder + "/null_pvalues.txt",
                    random_genesets(n_random_snps, slop),
                    min_observations,
                    max_iterations,
                    null_snpset_replicates
                );
            } else {
                // Calculate p-values for matched null gene sets.
                calculate_pvalues(
                    out_folder + "/null_pvalues.txt",
                    matched_genesets(),
                    min_observations,
                    max_iterations,
                    null_snpset_replicates
                );
            }
        }

        std::cout << timestamp() << " # done.\n"
                  << std::flush;
    }

    std::vector<std::vector<ulong> > genesets;
    for (auto item : _user_genesets) {
        genesets.push_back(item.second);
    }

    // Report p-values and gene identifiers for each SNP-column pair.
    report_pvalues(out_folder + "/snp_pvalues.txt", _user_genesets);

    std::cout << timestamp()
              << " # Computing one column at a time ...\n"
              << std::flush;

    // Calculate p-values for the user's SNP set.
    calculate_pvalues(
        out_folder + "/pvalues.txt",
        genesets,
        min_observations,
        max_iterations,
        1L
    );

    std::cout << timestamp() << " # done.\n\n";
}

void snpspec::write_args(
    std::string user_snpset_file,
    std::string gene_matrix_file,
    std::string gene_intervals_file,
    std::string snp_intervals_file,
    std::string null_snps_file,
    std::string condition_file,
    std::string out_folder,
    ulong slop,
    int threads,
    ulong null_snpset_replicates,
    ulong min_observations,
    ulong max_iterations,
    std::ostream & stream
)
{
    stream << "snpspec --snps " << user_snpset_file << "\n"
           << "        --gene-matrix " << gene_matrix_file << "\n"
           << "        --gene-intervals " << gene_intervals_file << "\n"
           << "        --snp-intervals " << snp_intervals_file << "\n"
           << "        --null-snps " << null_snps_file << "\n";

    if (condition_file.length() > 0) {
        stream << "        --condition " << condition_file << "\n";
    }

    stream << "        --out " << out_folder << "\n"
           << "        --slop " << slop << "\n"
           << "        --threads " << threads << "\n"
           << "        --null-snpsets " << null_snpset_replicates << "\n"
           << "        --min-observations " << min_observations << "\n"
           << "        --max-iterations " << max_iterations << "\n\n";
}

// Read an optionally gzipped text file and store the first column in a set of
// strings.
void snpspec::read_names(std::string filename, std::set<std::string> & names)
{
    gzifstream str(filename.c_str());
    if (!str.is_open()) {
        std::cerr << "ERROR: Cannot open " + filename << std::endl;
        exit(EXIT_FAILURE);
    }
    names.clear();
    Row row;
    while (str >> row) {
        names.insert(row[0]);
    }
    std::cout << timestamp() << " # \"" + filename + "\" has "
              << names.size() << " items." << std::endl;
}

// Given the name of a SNP, look up its interval and find overlapping genes.
// Report the offsets to lookup the genes in the gene matrix.
std::vector<ulong> snpspec::snp_geneset(std::string snp, ulong slop)
{
    auto snp_interval = _snp_intervals[snp];
    std::vector<ulong> indices;

    // Find overlapping genes.
    std::vector<Interval<ulong> > gene_intervals;
    _gene_interval_tree[snp_interval.chrom].findOverlapping(
        snp_interval.start,
        snp_interval.end,
        gene_intervals
    );

    if (gene_intervals.size() == 0) {
        _gene_interval_tree[snp_interval.chrom].findOverlapping(
            std::max(1UL, snp_interval.start - slop),
            snp_interval.end + slop,
            gene_intervals
        );
    }

    if (gene_intervals.size() > 0) {
        // Indices used for lookup in the gene matrix.
        for (auto interval : gene_intervals) {
            indices.push_back(interval.value);
        }
    }

    return indices;
}

// Generate a number of random SNPs given a filename like "random20".
void snpspec::random_snps(
    std::string filename,
    std::set<std::string> & names,
    ulong slop
)
{
    // Grab the desired length of the randomly generated SNP list.
    std::string::size_type sz;
    // Skip the first 6 characters of "random20".
    int n = std::stoi(filename.substr(6), &sz);

    // We can't access the elements of a set quickly, so just copy it. This
    // uses extra memory, but I have about 600K SNPs in this list for TGP
    // so it's not bad.
    static std::vector<std::string> null_snps = make_vector(_null_snp_names);

    // Standard Mersenne Twister random number generator.
    static std::mt19937 generator;

    // Clear out the old set of SNP names.
    names.clear();

    while (names.size() < n) {
        // Pick a random null SNP name.
        std::uniform_int_distribution<ulong>
        distribution(0, null_snps.size() - 1);

        auto r = distribution(generator);
        
        // Sanity check. The SNP name must be in our map.
        if (_snp_intervals.count(null_snps[r]) == 0) {
            continue;
        }

        // Grab the gene set for the given SNP.
        auto geneset = snp_geneset(null_snps[r], slop);

        // The SNP must overlap at least one gene.
        if (geneset.size() == 0) {
            continue;
        }

        // Put a new SNP name in.
        names.insert(null_snps[r]);
    }
}

// Read an optionally gzipped BED file and store the genomic intervals in
// a map of name => interval.
void snpspec::read_bed_intervals(
    std::string filename,
    std::map<std::string, genomic_interval> & intervals
)
{
    gzifstream stream(filename.c_str());
    if (!stream.is_open()) {
        std::cerr << "ERROR: Cannot open " + filename << std::endl;
        exit(EXIT_FAILURE);
    }
    BEDRow row;
    while (stream >> row) {
        intervals[row.name] = row.i;
    }
    std::cout << timestamp() << " # \"" + filename + "\" has "
              << intervals.size() << " intervals." << std::endl;
}

// Read an optionally gzipped BED file and store the genomic intervals in
// an interval tree. (Actually, one interval tree for each chromosome.)
void snpspec::read_bed_interval_tree(
    std::string filename,
    const std::vector<std::string> & row_names,
    std::map<std::string, IntervalTree<ulong> > & tree
)
{
    gzifstream stream(filename.c_str());
    if (!stream.is_open()) {
        std::cerr << "ERROR: Cannot open " + filename << std::endl;
        exit(EXIT_FAILURE);
    }

    // Convert the vector to a set.
    std::set<std::string> row_names_set(row_names.begin(), row_names.end());

    // Map a chromosome name to a vector of intervals.
    typedef Interval<ulong> interval;
    std::map<std::string, vector<interval> > intervals;

    // Rather than storing the gene identifiers in the tree, we'll store the
    // indices of the gene identifiers in the provided vector.
    std::map<std::string, ulong> index;
    for (ulong i = 0; i < row_names.size(); i++) {
        index[row_names.at(i)] = i;
    }

    ulong skipped_genes = 0;
    BEDRow row;
    while (stream >> row) {
        // Skip the gene if it is not present in the gene matrix.
        if (row_names_set.count(row.name) != 0) {
            // Add an interval to the vector for the corresponding chromosome.
            // (The value stored in the tree is a ulong that is an index to
            // the row names of the gene matrix. It is later retrieved with
            // the findOverlapping() method.)
            intervals[row.i.chrom].push_back(
                interval(row.i.start, row.i.end, index[row.name])
            );
        } else {
            skipped_genes++;
        }
    }

    std::cout << timestamp()
              << " # Skipped loading " << skipped_genes
              << " gene intervals because they are absent from the"
              << " --gene-matrix file."
              << std::endl;

    // Loop through the chromosomes.
    for (auto item : intervals) {
        // item.first is the name of a chromosome.
        tree[item.first] = IntervalTree<ulong> (intervals[item.first]);
    }
}

void snpspec::read_gct(
    std::string filename,
    std::vector<std::string> & row_names,
    std::vector<std::string> & col_names,
    MatrixXd & data
)
{
    gzifstream stream(filename.c_str());
    if (!stream.is_open()) {
        std::cerr << "ERROR: Cannot open " + filename << std::endl;
        exit(EXIT_FAILURE);
    }

    // Check that the first line is correct.
    std::string str;
    stream >> str;
    if (str.find("#1.2") != 0) {
        std::cerr << "ERROR: Not a GCT file " + filename << std::endl;
        exit(EXIT_FAILURE);
    }

    // Read the number of rows and columns.
    unsigned int rows, cols;
    stream >> rows >> cols;
    // Skip to next line.
    std::getline(stream, str);

    if (rows <= 0 || cols <= 0) {
        std::cerr
                << "ERROR: Line 2 of GCT file is malformed " + filename
                << std::endl;
        exit(EXIT_FAILURE);
    } else {
        std::cout
                << timestamp()
                << " # \"" + filename + "\" has "
                << rows << " rows, " << cols << " columns." << std::endl;
    }

    // Resize our matrix to hold all of the data.
    data.resize(rows, cols);

    // Skip "Name" and "Description".
    std::getline(stream, str, '\t');
    std::getline(stream, str, '\t');
    // Read the column names.
    for (int c = 0; c < cols - 1; c++) {
        std::getline(stream, str, '\t');
        str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
        col_names.push_back(str);
    }
    stream >> str;
    col_names.push_back(str);
    // Skip to next line.
    std::getline(stream, str);

    for (int r = 0; r < rows; r++) {
        // Read the Name in the first column.
        std::getline(stream, str, '\t');
        row_names.push_back(str);

        // Skip the Description column.
        std::getline(stream, str, '\t');

        // Read the data in this row.
        for (int c = 0; c < cols; c++) {
            stream >> data(r, c);
        }

        // Skip to next line.
        std::getline(stream, str);
    }
}

void snpspec::overlap_genes(
    std::set<std::string> & snp_names,
    std::set<std::string> & absent_snp_names,
    std::map<std::string, std::vector<ulong> > & genesets,
    std::vector<ulong> & geneset_sizes,
    ulong slop
)
{
    std::cout << timestamp()
              << " # Overlapping SNP intervals with gene intervals ...\n";

    // Clear out the old gene sets from previous runs.
    absent_snp_names.clear();
    geneset_sizes.clear();
    genesets.clear();

    for (auto snp : snp_names) {
        if (_snp_intervals.count(snp) == 0) {
            std::cout << timestamp() << " # "
                      << snp << " not found in --snp-intervals file.\n";
            absent_snp_names.insert(snp);
        } else {
            // Find the SNP's interval and find overlapping genes.
            std::vector<ulong> gene_ids = snp_geneset(snp, slop);
            if (gene_ids.size() > 0) {
                genesets[snp] = gene_ids;
            }
        }
    }
    std::cout << timestamp() << " # done.\n";
}

void snpspec::merge_user_snps(
    std::set<std::string> & snp_names,
    std::map<std::string, std::vector<ulong> > & genesets,
    std::vector<ulong> & geneset_sizes
)
{
    // Create new variables and fill them after merging SNPs.
    std::set<std::string> new_snp_names;
    std::map<std::string, std::vector<ulong> > new_genesets;
    std::vector<ulong> new_geneset_sizes;

    std::set<std::string> merged_snps;

    std::cout << timestamp()
              << " # Merging SNPs with shared genes ...\n";

    int count_snps = 0;
    int count_loci = 0;

    // Brute force, check all pairs of SNPs.
    // If two SNPs reside on the same chromosome and have shared or
    // overlapping genes, then merge them.
    for (auto a : snp_names) {
        if (genesets.count(a) == 0) continue;
        if (merged_snps.count(a) > 0) continue;

        std::vector<ulong> genes_a = genesets[a];
        std::sort(genes_a.begin(), genes_a.end());

        std::string merged_snp = a;

        for (auto b : snp_names) {
            if (a.compare(b) == 0) continue;
            if (genesets.count(b) == 0) continue;
            if (merged_snps.count(b) > 0) continue;

            std::vector<ulong> genes_b = genesets[b];
            std::sort(genes_b.begin(), genes_b.end());
            
            // Find the union of the two gene sets.
            std::vector<ulong> genes_ab (genes_a.size() + genes_b.size());
            std::vector<ulong>::iterator it;
            it = std::set_union(genes_a.begin(), genes_a.end(),
                                genes_b.begin(), genes_b.end(),
                                genes_ab.begin());
            genes_ab.resize(it - genes_ab.begin());

            if (genes_ab.size() < genes_a.size() + genes_b.size()) {
                if (count_snps == 0) {
                    count_snps = 2;
                } else {
                    count_snps += 2;
                }
                merged_snp += ":" + b;
                genes_a = genes_ab;
                merged_snps.insert(a);
                merged_snps.insert(b);
            }
        }

        if (merged_snp.find(":") != std::string::npos) {
            std::cout << timestamp()
                      << " # " << merged_snp << "\n";
            count_loci++;
        }

        new_snp_names.insert(merged_snp);
        new_genesets[merged_snp] = genes_a;
        new_geneset_sizes.push_back(genes_a.size());
    }

    snp_names = new_snp_names;
    genesets = new_genesets;
    geneset_sizes = new_geneset_sizes;

    if (count_snps > 0) {
        std::cout << timestamp() << " # done. Merged "
                  << count_snps << " SNPs into "
                  << count_loci << " loci.\n" << std::flush;
    } else {
        std::cout << timestamp() << " # done. No SNPs were merged.\n"
                  << std::flush;
    }
}

void snpspec::report_user_snp_genes(const std::string & filename)
{
    std::cout << timestamp() << " # Writing \"" + filename + "\" ...\n";

    ofstream stream(filename);

    // Print the column names.
    stream << "chrom\tstart\tend\tname\tn_genes\tgenes\n";

    // Print a row for each of the user's SNPs.
    for (auto snp : _user_snp_names) {
        // If the snp is not found in --snp-intervals then print NA.
        if (_user_genesets.count(snp) == 0) {
            stream << "NA\tNA\tNA\t" << snp << "\tNA\tNA\n";
        } else {
            std::string chrom;
            ulong start = 0;
            ulong end = 0;

            // This is a merged SNP, so print the interval that captures the
            // entire merged locus.
            if (snp.find(":") != std::string::npos) {
                // Loop through the snps and choose the min start, max end.
                for (auto merged_snp : split_string(snp, ':')) {
                    auto snp_interval = _snp_intervals[merged_snp];
                    chrom = snp_interval.chrom;
                    if (start == 0 || snp_interval.start < start) {
                        start = snp_interval.start;
                    }
                    if (end == 0 || snp_interval.end > end) {
                        end = snp_interval.end;
                    }
                }
            } else {
                auto snp_interval = _snp_intervals[snp];
                chrom = snp_interval.chrom;
                start = snp_interval.start;
                end = snp_interval.end;
            }

            auto geneset = _user_genesets[snp];
            // Print a BED line with two extra columns:
            //      - number of overlapping genes
            //      - Entrez IDs
            stream << chrom << '\t'
                   << start << '\t'
                   << end << '\t'
                   << snp << '\t'
                   << geneset.size() << '\t';

            if (geneset.size() > 0) {
                // Print the first gene, then prepend a comma to the next.
                stream << _row_names.at(geneset.at(0));

                for (int i = 1; i < geneset.size(); i++) {
                    stream << ',' << _row_names.at(geneset.at(i));
                }
            }
            stream << std::endl;
        }
    }
    stream.close();

    std::cout << timestamp() << " # done." << std::endl;
}

/*
// Drop intervals for SNPs absent from "--null-snps". This is desirable to
// reduce memory usage, but must be done after the user's SNPs are assigned
// intervals.
void snpspec::drop_snp_intervals()
{
    ulong dropped_snps = 0;
    auto it = _snp_intervals.begin();
    while (it != _snp_intervals.end()) {
        if (_null_snp_names.count(it->first) == 0) {
            _snp_intervals.erase(it++);
            dropped_snps++;
        } else {
            ++it;
        }
    }
    std::cout << timestamp()
              << " # Dropped " << dropped_snps
              << " SNP intervals that are absent from the provided null set."
              << std::endl;
}
*/

void snpspec::report_missing_conditions()
{
    if (_condition_names.size() == 0) {
        return;
    }
    std::set<std::string> _col_names_set(_col_names.begin(), _col_names.end());
    std::set<std::string> _condition_difference;

    std::set_difference(
        _condition_names.begin(), _condition_names.end(),
        _col_names_set.begin(), _col_names_set.end(),
        std::inserter(_condition_difference, _condition_difference.begin())
    );

    if (_condition_difference.size() > 0) {
        std::cerr << "ERROR: Conditions not found in --gene-matrix file:"
                  << std::endl;
        for (auto name : _condition_difference) {
            std::cerr << name << std::endl;
        }
        exit(EXIT_FAILURE);
    }
}

// Condition the gene matrix on the specified column names. Each column is
// projected onto a condition column, and its projection is substracted.
void snpspec::condition(
    MatrixXd & matrix,
    std::set<std::string> & col_names
)
{
    std::vector<std::string>
    new_col_names (_col_names.begin(), _col_names.end());

    std::vector<size_t> idxs;
   
    for (const auto & col_name : col_names) {
        auto it = std::find(_col_names.begin(), _col_names.end(), col_name);
        size_t col_index = it - _col_names.begin();

        idxs.push_back(col_index);

        const auto b = VectorXd(matrix.col(col_index));

        for (size_t col = 0; col < _col_names.size(); col++) {
            auto a = matrix.col(col);
            auto projection = a.dot(b) / b.dot(b) * b;
            
            matrix.col(col) -= projection;
        }
    }
    // Remove the condition columns from the matrix.
    removeColumns(idxs, matrix);
    // Remove the condition column names. Sort descending, then delete them.
    std::sort(idxs.begin(), idxs.end(), std::greater<int>());
    for (auto idx : idxs) {
        new_col_names.erase(new_col_names.begin() + idx);
    }
    // Replace the old column names with new ones.
    _col_names = new_col_names;
}

void snpspec::bin_genesets(ulong slop, ulong max_genes)
{
    for (const auto & item : _snp_intervals) {
        std::string snp = item.first;

        // We want to sample from the list in "--null-snps".
        if (_null_snp_names.count(snp) == 0) continue;

        std::vector<ulong> geneset = snp_geneset(snp, slop);

        // Put the geneset in a bin that corresponds to its size.
        ulong n_genes = geneset.size();
        if (n_genes > 0) {
            // Put an upper limit on the number of genes in a set. So, if
            // a geneset actually has more genes, that's ok.
            if (n_genes > max_genes) {
                n_genes = max_genes;
            }
            // Indices used for lookup in the gene matrix.
            _geneset_bins[n_genes].push_back(geneset);
        }
    }
}

// Generate a vector of vectors. Each inner vector contains gene indices for
// looking up rows in the gene matrix.
std::vector<std::vector<ulong> > snpspec::matched_genesets()
{
    // Standard Mersenne Twister random number generator.
    static std::mt19937 generator;

    std::vector<std::vector<ulong> > genesets;
    for (auto s : _user_geneset_sizes) {
        // Uniform integer distribution.
        std::uniform_int_distribution<ulong>
        distribution(0, _geneset_bins[s].size() - 1);

        auto r = distribution(generator);

        genesets.push_back(_geneset_bins[s].at(r));
    }
    return genesets;
}

// Same as matched_genesets(), but pick gene sets randomly without matching.
std::vector<std::vector<ulong> > snpspec::random_genesets(int n, ulong slop)
{
    std::vector<std::vector<ulong> > genesets;
    // Get a list of random SNPs.
    std::set<std::string> snps;
    random_snps("random" + std::to_string(n), snps, slop);

    for (auto snp : snps) {
        auto geneset = snp_geneset(snp, slop);
        if (geneset.size() > 0) {
            genesets.push_back(geneset);
        }
    }
    return genesets;
}

// Returns a score for a column in the binary gene matrix using the
// genes in the given set of gene sets.
double snpspec::score_binary(
    const ulong & col,
    const std::vector<std::vector<ulong> > & genesets
)
{
    ulong n = _binary_sums(col);
    double p = _binary_probs(col);
    double score = 0.0;
    for (auto geneset : genesets) {
        int k = 0;
        for (auto gene_id : geneset) {
            if (_gene_matrix(gene_id, col) > 0) {
                k++;
            }
        }
        score += -log10(gsl_ran_binomial_pdf(k, p, n));
    }
    return std::isfinite(score) ? score : 0.0;
}

// Returns a score for a column in the quantitative gene matrix using
// the genes in the given set of gene sets.
double snpspec::score_quantitative(
    const ulong & col,
    const std::vector<std::vector<ulong> > & genesets
)
{
    double score = 0.0;
    for (auto geneset : genesets) {
        double percentile = 1.0;
        for (auto gene_id : geneset) {
            percentile = std::min(percentile, _gene_matrix(gene_id, col));
        }
        if (percentile < 1.0) {
            // Each gene set contributes to the score.
            score += -log10(1 - pow(1 - percentile, geneset.size()));
        }
    }
    return std::isfinite(score) ? score : 0.0;
}


void snpspec::report_pvalues(
    const std::string filename,
    const std::map<std::string, std::vector<ulong> > genesets
)
{
    std::cout << timestamp() << " # Writing \"" + filename + "\" ...\n";

    // Open the file.
    ofstream stream (filename);
    // Write the header.
    stream << "marker\tcolumn\tgene\tpvalue\n";
    // Iterate through each SNP's gene set.
    for (const auto & kv : genesets) {
        // Iterate through each column.
        for (int col = 0; col < _col_names.size(); col++) {
            double pvalue = 1;
            std::string min_gene;
            if (_binary_gene_matrix) {
                ulong n = _binary_sums(col);
                double p = _binary_probs(col);
                int k = 0;
                for (auto gene_id : kv.second) {
                    if (_gene_matrix(gene_id, col) > 0) {
                        k++;
                    }
                }
                pvalue = gsl_ran_binomial_pdf(k, p, n);
                // The p-value does not depend on a single gene for binary
                // matrices, but instead for the whole gene set.
                min_gene = "NULL";
            } else {
                double percentile = 1.0;
                for (auto gene_id : kv.second) {
                    if (_gene_matrix(gene_id, col) < percentile) {
                        percentile = _gene_matrix(gene_id, col);
                        min_gene = _row_names[gene_id];
                    }
                }
                if (percentile < 1.0) {
                    // Each gene set contributes to the score.
                    pvalue = 1 - pow(1 - percentile, kv.second.size());
                }
            }
            // The SNP's name, column name, best rank gene, p-value.
            stream << kv.first << "\t"
                   << _col_names[col] << "\t"
                   << min_gene << "\t"
                   << pvalue << "\n";
        }
    }
    stream.close();

    std::cout << timestamp() << " # done." << std::endl;
}

void snpspec::calculate_pvalues(
    std::string filename,
    std::vector<std::vector<ulong> > genesets,
    long min_observations,
    long max_iterations,
    long replicates
)
{
    static int replicate = -1;
    replicate++;

    // Set the appropriate scoring function.
    auto score_function = &snpspec::score_quantitative;
    if (_binary_gene_matrix) {
        score_function = &snpspec::score_binary;
    }

    std::fstream stream;
    if (replicates <= 1) {
        // Overwrite the file.
        stream.open(filename, std::fstream::out);
        stream << "name\tpvalue\tnulls_observed\tnulls_tested" << std::endl;
    } else {
        // Append to the file.
        stream.open(filename, std::fstream::out | std::fstream::app);
    }

    for (ulong col = 0; col < _gene_matrix.cols(); col++) {
        // Shared across all threads.
        double user_score = (*this.*score_function)(col, genesets);

        // The user's SNPs scored 0, so don't bother testing.
        if (user_score <= 0) {
            stream << _col_names.at(col) << "\t1.0\t0\t0";
            if (replicates > 1) {
                stream << "\t" << replicate;
            }
            stream << std::endl;
            continue;
        }

        long nulls_tested = 0;
        long nulls_observed = 0;

        for (auto count : iterations(100, max_iterations)) {
            #pragma omp parallel
            {
                // Private to each thread.
                long thread_observed = 0;

                // Each thread will complete some fraction of this loop.
                #pragma omp for
                for (long i = 0; i < count; i++) {
                    // Call the appropriate scoring function.
                    if ((*this.*score_function)(col, matched_genesets())
                            >= user_score) {
                        thread_observed += 1;
                    }
                }

                // Each thread counts its own results, then we sum them.
                #pragma omp critical
                {
                    nulls_observed += thread_observed;
                }
            }
            // Count how many total iterations we performed.
            nulls_tested += count;

            // A null SNP set scored higher than the user's SNP set enough
            // times that we are confident in the column's p-value.
            if (nulls_observed >= min_observations) {
                break;
            }
        }

        double pvalue = double(nulls_observed) / double(nulls_tested);

        if (nulls_observed == 0) pvalue = 1.0 / double(nulls_tested);

        stream << _col_names.at(col) << '\t' << pvalue << '\t'
               << nulls_observed << '\t' << nulls_tested;

        // Display a period for each column.
        if (replicates > 1) {
            stream << '\t' << replicate;
        } else {
            std::cout << '.';
            if ((col + 1) %  5 == 0) std::cout << ' ';
            if ((col + 1) % 10 == 0) std::cout << ' ';
            if ((col + 1) % 50 == 0) std::cout << col + 1 << '\n';
            std::cout << std::flush;
        }

        stream << '\n' << std::flush;
    }

    // Display a period for each replicate.
    if (replicates > 1) {
        std::cout << '.';
        if ((replicate + 1) %  5 == 0) std::cout << ' ';
        if ((replicate + 1) % 10 == 0) std::cout << ' ';
        if ((replicate + 1) % 50 == 0) std::cout << replicate + 1 << '\n';
        std::cout << std::flush;
    } else {
        std::cout << '\n' << std::flush;
    }

    stream.close();
}
