#include "snpspec.h"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 0
#define omp_set_num_threads() 0
#endif

snpspec::snpspec(
    std::vector<std::string> user_snpset_files,
    std::string expression_file,
    std::string gene_intervals_file,
    std::string snp_intervals_file,
    std::string null_snps_file,
    std::string condition_file,
    std::string out_folder,
    int slop,
    int threads,
    long null_snpset_replicates,
    long min_observations,
    long max_iterations 
)
{
    std::cout
            << "snpspec --snps\n";

    for (int i = 0; i < user_snpset_files.size() - 1; i++) {
        std::cout << "              " << user_snpset_files.at(i) << ",\n";
    }
    std::cout << "              " << user_snpset_files.back() << "\n";

    std::cout
            << "        --expression " + expression_file + "\n"
            << "        --gene-intervals " + gene_intervals_file + "\n"
            << "        --snp-intervals " + snp_intervals_file + "\n"
            << "        --null-snps " + null_snps_file + "\n";

    if (condition_file.length() > 0) {
        std::cout << "        --condition " + condition_file + "\n";
    }

    std::cout
            << "        --out " + out_folder + "\n"
            << "        --slop " << slop << "\n"
            << "        --threads " << threads << "\n"
            << "        --null-snpsets " << null_snpset_replicates << "\n"
            << "        --min-observations " << min_observations << "\n"
            << "        --max-iterations " << max_iterations << "\n"
            << std::endl;

    // Read names.
    std::cout << timestamp() << " # Reading files..." << std::endl;
    read_names(null_snps_file, _null_snp_names);

    // Optional.
    if (condition_file.length() > 0) {
        read_names(condition_file, _condition_names);
    }

    // Read SNP names and intervals.
    read_bed_intervals(snp_intervals_file, _snp_intervals);

    // Read the gene expression GCT file.
    read_gct(expression_file, _row_names, _col_names, _expression);

    // Read the gene intervals but only keep the ones listed in the GCT.
    read_bed_interval_tree(
        gene_intervals_file,
        _row_names,
        _gene_interval_tree
    );

    std::cout << timestamp() << " # done." << std::endl;

    // Report names from the conditions file that are absent from the
    // expression file.
    report_missing_conditions();

    // Drop all SNP intervals except those in the null set.
    //drop_snp_intervals();

    // Check if the matrix is binary by reading the first column.
    if (is_binary(_expression.col(0))) {
        // Let the user know we detected it.
        std::cout << timestamp() << " # Expression is binary." << std::endl;
        // Cache these values ahead of time.
        _binary_sums = _expression.colwise().sum();
        _binary_probs = _binary_sums / _expression.rows();
    } else {
        // TODO Condition the matrix on the specified columns.

        // Normalize the matrix.
        _expression.colwise().normalize();

        // Reverse percentile rank each column of the matrix.
        // So, a small value like 0.02 means the given gene is highly specific
        // to the column. A large value means the gene is non-specific.
        for (int i = 0; i < _expression.cols(); i++) {
            _expression.col(i) =
                rankdata(_expression.col(i)) / _expression.rows();
        }
    }

    // 1. Find a geneset for each SNP by querying the gene interval tree.
    // 2. Bin genesets by size. (This will be used to generate SNP sets.)
    const int MAX_GENES = 10;
    bin_genesets(slop, MAX_GENES);

    // Check for enrichment of each column in parallel.
    // Ensure that a valid number of threads is used.
    threads = clamp(threads, 1, cpu_count());
    omp_set_num_threads(threads);

    for (auto f : user_snpset_files) {
        read_names(f, _user_snp_names);

        std::string path = out_folder + "/" + timestamp("%F_%H-%M-%S");
        mkpath(path);

        // TODO Merge SNPs that share genes or have overlapping genes.
//        for (int i = 0; i < _user_snp_names.size(); i++) {
//            for (int j = i + 1; j < _user_snp_names.size(); j++) {
//                
//            }
//        }

        // 1. Report the genes overlapping the user's SNPs.
        // 2. Clear and fill _user_genesets and _user_geneset_sizes.
        report_user_snp_genes(path + "/snp_genes.txt", slop);

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
                std::cout
                          << " from a pool of size " << item.second.size()
                          << std::endl;
            }
        }

        std::cout << timestamp()
                  << " # Computing up to " << scientific << max_iterations
                  << " iterations for each column in the expression with "
                  << fixed << threads << " threads...\n";
        std::cout << std::flush;

        // TODO This is an easter egg. Make an actual option.
        if (min_observations == 27) {
            calculate_pvalues(
                path + "/pvalues_null.txt",
                generate_snpset(),
                min_observations,
                max_iterations,
                null_snpset_replicates
            );
        } else {
            calculate_pvalues(
                path + "/pvalues.txt",
                _user_genesets,
                min_observations,
                max_iterations,
                1
            );
        }

        std::cout << timestamp() << " # done.\n\n";
    }
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
              << intervals.size() << " items." << std::endl;
}

// Read an optionally gzipped BED file and store the genomic intervals in
// an interval tree. (Actually, one interval tree for each chromosome.)
void snpspec::read_bed_interval_tree(
    std::string filename,
    const std::vector<std::string> & whitelist,
    std::map<std::string, IntervalTree<size_t> > & tree
)
{
    gzifstream stream(filename.c_str());
    if (!stream.is_open()) {
        std::cerr << "ERROR: Cannot open " + filename << std::endl;
        exit(EXIT_FAILURE);
    }

    // Convert the vector to a set.
    std::set<std::string> whiteset(whitelist.begin(), whitelist.end());

    // Map a chromosome name to a vector of intervals.
    typedef Interval<size_t> interval;
    std::map<std::string, vector<interval> > intervals;

    // Rather than storing the Entrez IDs in the tree, we'll store the indices
    // of the Entrez IDs in the provided vector.
    std::map<std::string, size_t> index;
    for (int i = 0; i < whitelist.size(); i++) {
        index[whitelist.at(i)] = i;
    }

    int skipped_genes = 0;
    BEDRow row;
    while (stream >> row) {
        // Add an interval to the vector for the corresponding chromosome.
        // (The value stored in the tree is a pair of integers, and it is
        // later retrieved with the findOverlapping() method.)
        if (whiteset.count(row.name) != 0) {
            intervals[row.i.chrom].push_back(
                interval(row.i.start, row.i.end, index[row.name])
            );
        } else {
            skipped_genes++;
        }
    }

    std::cout
            << timestamp()
            << " # Skipped loading " << skipped_genes
            << " gene intervals because they are absent from the"
            << " expression file."
            << std::endl;

    // Loop through the chromosomes.
    for (auto item : intervals) {
        // item.first is the name of a chromosome.
        tree[item.first] = IntervalTree<size_t> (intervals[item.first]);
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
    for (int c = 0; c < 2; c++) {
        stream >> str;
    }
    // Read the column names.
    for (int c = 0; c < cols; c++) {
        stream >> str;
        col_names.push_back(str);
    }

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

void snpspec::report_user_snp_genes(const std::string & filename, int slop)
{
    ofstream stream(filename);

    std::cout << timestamp() << " # Writing \"" + filename + "\" ...\n";

    // Print the column names.
    stream << "chrom\tstart\tend\tname\tn_genes\tgenes\n";

    // Clear out the old gene sets from previous runs.
    _user_geneset_sizes.clear();
    _user_genesets.clear();
    // TODO This would be nice to report to the user...
    //_user_genesets_absent.clear();

    // Print a row for each of the user's SNPs.
    for (auto snp : _user_snp_names) {
        // If snp's not found in the reference SNP intervals, print a blank.
        if (_snp_intervals.count(snp) == 0) {
            stream << "NA\tNA\tNA\t" << snp << "\tNA\tNA\n";
        } else {
            // Find its interval and find overlapping genes.
            genomic_interval snp_interval = _snp_intervals[snp];
            std::vector<Interval<size_t> > gene_intervals;
            _gene_interval_tree[snp_interval.chrom].findOverlapping(
                snp_interval.start,
                snp_interval.end,
                gene_intervals
            );

            // Adjust with slop and try again if we found no genes.
            if (gene_intervals.size() == 0) {
                _gene_interval_tree[snp_interval.chrom].findOverlapping(
                    std::max(1, snp_interval.start - slop),
                    snp_interval.end + slop,
                    gene_intervals
                );
            }

            // Print a BED line with two extra columns:
            //      - number of overlapping genes
            //      - Entrez IDs
            stream << snp_interval.chrom << '\t'
                   << snp_interval.start << '\t'
                   << snp_interval.end << '\t'
                   << snp << '\t'
                   << gene_intervals.size() << '\t';

            if (gene_intervals.size() > 0) {
                // Indices for lookup in the expression matrix.
                std::vector<size_t> gene_ids;

                // Save the sizes of the user's SNP's gene set. These values
                // are later used to generate random matched SNP sets.
                _user_geneset_sizes.push_back(gene_intervals.size());

                // Print the first gene, then prepend a comma to the next.
                stream << _row_names.at(gene_intervals.at(0).value);
                gene_ids.push_back(gene_intervals.at(0).value);

                for (int i = 1; i < gene_intervals.size(); i++) {
                    stream << ',' << _row_names.at(gene_intervals.at(i).value);
                    gene_ids.push_back(gene_intervals.at(i).value);
                }

                _user_genesets.push_back(gene_ids);
            }
            stream << std::endl;
        }
    }
    stream.close();

    std::cout << timestamp() << " # done.\n";
}

void snpspec::drop_snp_intervals()
{
    int dropped_snps = 0;
    auto it = _snp_intervals.begin();
    while (it != _snp_intervals.end()) {
        if (_null_snp_names.count(it->first) == 0) {
            _snp_intervals.erase(it++);
            dropped_snps++;
        } else {
            ++it;
        }
    }
    std::cout
            << timestamp()
            << " # Dropped " << dropped_snps
            << " SNP intervals that do not belong to the provided null set."
            << std::endl;
}

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
        std::cerr << "ERROR: Conditions not found in expression file:"
                  << std::endl;
        for (auto name : _condition_difference) {
            std::cerr << name << std::endl;
        }
        exit(EXIT_FAILURE);
    }
}

void snpspec::bin_genesets(int slop, int max_genes)
{
    //auto geneset_sizes = make_set(_user_geneset_sizes);
    for (auto item : _snp_intervals) {
        // Find overlapping genes.
        std::vector<Interval<size_t> > genes;
        _gene_interval_tree[item.second.chrom].findOverlapping(
            item.second.start,
            item.second.end,
            genes
        );

        if (genes.size() == 0) {
            _gene_interval_tree[item.second.chrom].findOverlapping(
                std::max(1, item.second.start - slop),
                item.second.end + slop,
                genes
            );
        }

        // Put the geneset in a bin that corresponds to its size.
        int n_genes = genes.size();
        if (n_genes > 0) {
            // Put an upper limit on the number of genes in a set. So, if
            // a geneset actually has more genes, that's ok.
            if (n_genes > max_genes) {
                n_genes = max_genes;
            }
            // We only care to maintain genesets that have the same number of
            // genes as the genesets that correspond to the user's SNPs.
            //if (geneset_sizes.count(n_genes) == 0) {
            //    continue;
            //}
            // Indices used for lookup in the expression.
            std::vector<size_t> indices;
            for (auto interval : genes) {
                indices.push_back(interval.value);
            }
            _geneset_bins[n_genes].push_back(indices);
        }
    }
}

// Generate a vector of vectors. Each inner vector contains gene indices for
// looking up rows in the expression.
std::vector<std::vector<size_t> > snpspec::generate_snpset()
{
    std::vector<std::vector<size_t> > snpset;
    for (auto s : _user_geneset_sizes) {
        int r = std::rand() % _geneset_bins[s].size();
        snpset.push_back(_geneset_bins[s].at(r));
    }
    return snpset;
}

MatrixXd snpspec::geneset_pvalues_binary(std::vector<size_t> & geneset)
{
    //MatrixXd m(_user_geneset_sizes.size(), _expression.cols());
    MatrixXd m(1, _expression.cols());

    // TODO Write the body of this function. Everything compiles and works!
    m(0, 0) = gsl_ran_binomial_pdf(3, 0.001, 20);

    return m;
}

// Returns a score for a column in the binary _expression matrix using the
// genes in the given set of gene sets.
double snpspec::score_binary(
    const size_t & col,
    const std::vector<std::vector<size_t> > & snpset
)
{
    int n = (int) _binary_sums(col);
    double p = _binary_probs(col);
    double score = 0.0;
    for (auto geneset : snpset) {
        int k = 0;
        for (auto gene_id : geneset) {
            if (_expression(gene_id, col) > 0) {
                k++;
            }
        }
        score += -log10(gsl_ran_binomial_pdf(k, p, n));
    }
    return std::isfinite(score) ? score : 0.0;
}

// Returns a score for a column in the quantitative _expression matrix using
// the genes in the given set of gene sets.
double snpspec::score_quantitative(
    const size_t & col,
    const std::vector<std::vector<size_t> > & snpset
)
{
    double score = 0.0;
    for (auto geneset : snpset) {
        double percentile = 1.0;
        for (auto gene_id : geneset) {
            percentile = std::min(percentile, _expression(gene_id, col));
        }
        if (percentile < 1.0) {
            // Each gene set contributes to the score.
            score += -log10(1 - pow(1 - percentile, geneset.size()));
        }
    }
    return std::isfinite(score) ? score : 0.0;
}


void snpspec::calculate_pvalues(
    std::string filename,
    std::vector<std::vector<size_t> > genesets,
    long min_observations,
    long max_iterations,
    long replicates
)
{
    // Set the appropriate scoring function.
    auto score_function = &snpspec::score_quantitative;
    if (_binary_expression) {
        score_function = &snpspec::score_binary;
    }

    ofstream stream(filename);
    stream << "name\tpvalue\tnulls_observed\tnulls_tested\n";
    stream << std::flush;

    for (int col = 0; col < _expression.cols(); col++) {
        // Shared across all threads.
        double user_score = (*this.*score_function)(col, genesets);

        // The user's SNPs scored 0, so don't bother testing.
        if (user_score <= 0) {
            stream << _col_names.at(col) << "\t1.0\t0\t0\n";
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
                    if ((*this.*score_function)(col, generate_snpset())
                            >= user_score) {
                        thread_observed += 1;
                    }
                }

                // Each thread counts its own results and we sum afterwards.
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

        double pvalue = (double) nulls_observed / (double) nulls_tested;

        stream << _col_names.at(col) << '\t' << pvalue << '\t'
               << nulls_observed << '\t' << nulls_tested << '\n';
        stream << std::flush;
    }

    stream.close();
}
