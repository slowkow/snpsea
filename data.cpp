#include "snpspec.h"

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 0
    #define omp_set_num_threads() 0
#endif

snpspec::snpspec(
    std::string user_snps_file,
    std::string expression_file,
    std::string gene_intervals_file,
    std::string snp_intervals_file,
    std::string null_snps_file,
    std::string condition_file,
    std::string out_folder,
    int slop,
    int processes,
    int min_observations,
    int permutations
) {
    std::cout
        << timestamp() << " # Started with arguments:\n"
        << "snpspec --snps " + user_snps_file + " \n"
        << "        --expression " + expression_file + "\n"
        << "        --gene-intervals " + gene_intervals_file + "\n"
        << "        --snp-intervals " + snp_intervals_file + "\n"
        << "        --null-snps " + null_snps_file + "\n"
        << "        --condition " + condition_file + "\n"
        << "        --out " + out_folder + "\n"
        << "        --slop " << slop << "\n"
        << "        --processes " << processes << "\n"
        << "        --min-observations " << min_observations << "\n"
        << "        --permutations " << permutations << "\n"
        << std::endl;

    // Read names.
    std::cout << timestamp() << " # Reading files..." << std::endl;
    read_names(user_snps_file, _user_snp_names);
    read_names(null_snps_file, _null_snp_names);
    read_names(condition_file, _condition_names);

    // Read SNP names and intervals.
    read_bed_intervals(snp_intervals_file, _snp_intervals);

    

    // Read the gene expression GCT file.
    read_gct(expression_file, _row_names, _col_names, _expression);
    std::cout << timestamp() << " # done." << std::endl;

    // Read the gene intervals but only keep the ones listed in the GCT.
    read_bed_interval_tree(
            gene_intervals_file,
            _row_names,
            _gene_interval_tree
    );

    // Report the genes overlapping the user's SNPs.
    report_user_snp_genes(out_folder + "/snp_genes.txt");
    
    // Drop all SNP intervals except those in the null set.
    drop_snp_intervals();

    // Report names from the conditions file that are absent from the
    // expression file.
    report_missing_conditions();

    // Check if the matrix is binary by reading the first column.
    if (is_binary(_expression.col(0))) {
        std::cout << timestamp() << " # Expression is binary." << std::endl;
        _binary_sums = _expression.colwise().sum();
        _binary_probs = _binary_sums / _expression.rows();
    } else {
        // Condition the matrix on the specified columns.
        
        // Normalize the matrix.
        _expression.colwise().normalize();

        // Percentile rank each column of the matrix.
        for (int i = 0; i < _expression.cols(); i++) {
            _expression.col(i) = rankdata(_expression.col(i));
        }
    }

    // Find a geneset for each SNP by querying the gene interval tree.
    // Bin genesets by size. This will be used to generate SNP sets later.
    bin_genesets();

    // Delete large variables that are no longer needed.
    //delete[] & _snp_intervals;
    //delete[] & _gene_interval_tree;

    // Show me some random SNP sets.
    /*
    for (int i = 0; i < 10; i++) {
        auto snpset = generate_snpset();
        for (auto it = snpset.begin(); it != snpset.end(); it++) {
            std::cout << it->size() << " { ";
            for (int j = 0; j < it->size(); j++ ) {
                std::cout << it->at(j) << " ";
            }
            std::cout << "}" << std::endl;
        }
        std::cout << std::endl;
    }
    */

    // Check for enrichment of each column in parallel.
    // Ensure that a valid number of processes is used.
    processes = clamp(processes, 1, cpu_count());
    omp_set_num_threads(processes);

    std::cout << timestamp() 
        << " # Computing scores for null SNP sets with " 
        << processes << " threads...\n";
    std::cout << std::flush;

    ofstream stream(out_folder + "/pvalues.txt");
    stream << "name\tpvalue\tnulls_observed\tnulls_tested\n";

    for (int i = 0; i < _expression.cols(); i++) {
        // Shared across all threads.
        double user_score = score_binary(i, _user_genesets);

        // The user's SNPs scored 0, so don't bother testing.
        if (user_score <= 0) {
            stream << _col_names.at(i) << "\t1.0\t0\t0\n";
            continue;
        }

        //std::vector<int> nulls_observed(processes);
        int nulls_tested = 0;
        int nulls_observed = 0;

        const int LOOPS = 1000;
        #pragma omp parallel
        {
            // Private to each thread.
            //int thread_id = omp_get_thread_num();
            int thread_observed = 0;

            // Each thread counts its own results and we sum them afterwards.
            #pragma omp for
            for (int j = 0; j < LOOPS; j++) {
                if (score_binary(i, generate_snpset()) >= user_score) {
                    //nulls_observed.at(thread_id) += 1;
                    thread_observed += 1;
                }
            }

            #pragma omp critical
            {
                nulls_observed += thread_observed;
            }
        }

        nulls_tested += LOOPS;
        //int observed = 0;
        //for (int j = 0; j < processes; j++) {
        //    observed += nulls_observed.at(j);
        //}
        double pvalue = (double) nulls_observed / (double) nulls_tested;

        stream << _col_names.at(i) << '\t' << pvalue << '\t'
            << nulls_observed << '\t' << nulls_tested << '\n';
    }

    stream.close();

    std::cout << timestamp() << " # done." << std::endl;
    //MatrixXd geneset_pvalues = geneset_pvalues_binary(generate_snpset().at(0));
}

// Read an optionally gzipped text file and store the first column in a set of
// strings.
void snpspec::read_names(std::string filename, std::set<std::string> & names)
{
    gzifstream str(filename.c_str());
    if(!str.is_open()) {
        std::cerr << "ERROR: Cannot open " + filename << std::endl;
        exit(EXIT_FAILURE);
    }
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
) {
    gzifstream stream(filename.c_str());
    if(!stream.is_open()) {
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
) {
    gzifstream stream(filename.c_str());
    if(!stream.is_open()) {
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
        << " gene intervals because they are absent from the expression file."
        << std::endl;

    // Loop through the chromosomes.
    for (auto it = intervals.begin(); it != intervals.end(); it++) {
        // it->first is the name of a chromosome.
        tree[it->first] = IntervalTree<size_t>(intervals[it->first]);
    }
}

void snpspec::read_gct(
        std::string filename,
        std::vector<std::string> & row_names,
        std::vector<std::string> & col_names,
        MatrixXd & data
) {
    gzifstream stream(filename.c_str());
    if(!stream.is_open()) {
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

void snpspec::report_user_snp_genes(const std::string & filename)
{
    ofstream stream(filename);
    // Print a row for each of the user's SNPs.
    for (auto it = _user_snp_names.begin();
            it != _user_snp_names.end(); it++) {
        // If it's not found in the reference SNP intervals, print a blank.
        if (_snp_intervals.count(*it) == 0) {
            stream << "NA\tNA\tNA\t" << *it << "\tNA\tNA\n";
        } else {
            // Find its interval and find overlapping genes.
            genomic_interval i = _snp_intervals[*it];
            std::vector<Interval<size_t> > genes;
            _gene_interval_tree[i.chrom].findOverlapping(i.start, i.end, genes);
            // Print a BED line with one two extra columns:
            //      - number of overlapping genes
            //      - Entrez IDs
            stream
                << i.chrom << '\t' << i.start << '\t' << i.end << '\t'
                << *it << '\t' << genes.size() << '\t';
            if (genes.size() > 0) {
                // Save the sizes of the user's SNP's gene set. These values
                // are later used to generate random matched SNP sets.
                _user_snp_geneset_sizes.push_back(genes.size());
                std::vector<size_t> gene_ids;
                // Print the first gene, then prepend a comma to the next.
                stream << _row_names.at(genes.at(0).value);
                for (auto it2 = genes.begin() + 1; it2 != genes.end(); it2++) {
                    stream << ',' << _row_names.at(it2->value);
                    gene_ids.push_back(it2->value);
                }
                _user_genesets.push_back(gene_ids);
            }
            stream << std::endl;
        }
    }
    stream.close();
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
        for (auto it = _condition_difference.begin();
                it != _condition_difference.end(); it++) {
            std::cerr << *it << std::endl;
        }
        exit(EXIT_FAILURE);
    }
}

void snpspec::bin_genesets() {
    const int MAX_GENES = 10;
    for (int i = 0; i < _user_snp_geneset_sizes.size(); i++) {
        if (_user_snp_geneset_sizes.at(i) > MAX_GENES) {
            _user_snp_geneset_sizes.at(i) = MAX_GENES;
        }
    }
    auto geneset_sizes = make_set(_user_snp_geneset_sizes);
    for (auto it = _snp_intervals.begin(); it != _snp_intervals.end(); it++) {
        // Find overlapping genes.
        std::vector<Interval<size_t> > genes;
        _gene_interval_tree[it->second.chrom].findOverlapping(
                it->second.start, it->second.end, genes
        );
        int n_genes = genes.size();
        if (n_genes > 0) {
            // Put an upper limit on the number of genes in a set. So, if
            // a genset actually has more genes, that's ok.
            if (n_genes > MAX_GENES) {
                n_genes = MAX_GENES;
            }
            // We only care to maintain genesets that have the same number of
            // genes as the genesets that correspond to the user's SNPs.
            if (geneset_sizes.count(n_genes) == 0) {
                continue;
            }
            std::vector<size_t> indices;
            for (int i = 0; i < genes.size(); i++) {
                indices.push_back(genes[i].value);
            }
            _geneset_bins[n_genes].push_back(indices);
        }
    }
    // Report how many genesets exist of each size.
    for (auto it = _geneset_bins.begin(); it != _geneset_bins.end(); it++) {
        std::cout << timestamp()
            << " # Gene sets with size " << it->first << ": "
            << it->second.size() << std::endl;
    }
}

// Generate a vector of vectors. Each inner vector contains gene indices for
// looking up rows in the expression.
std::vector<std::vector<size_t> > snpspec::generate_snpset() {
    std::vector<std::vector<size_t> > snpset;
    for (int i = 0; i < _user_snp_geneset_sizes.size(); i++) {
        int s = _user_snp_geneset_sizes.at(i);
        int r = std::rand() % _geneset_bins[s].size();
        snpset.push_back(_geneset_bins[s].at(r));
    }
    return snpset;
}

MatrixXd snpspec::geneset_pvalues_binary(std::vector<size_t> & geneset) {
    //MatrixXd m(_user_snp_geneset_sizes.size(), _expression.cols());
    MatrixXd m(1, _expression.cols());

    // TODO Write the body of this function. Everything compiles and works!
    m(0, 0) = gsl_ran_binomial_pdf(3, 0.001, 20);

    return m;
}

// One function. This way we progress one column at a time.
//      Generate a snpset.
//      Calculate a score for a single column.
//      Compare the score to the reference.
//      Return 1 or 0.
double snpspec::score_binary(
        const size_t & col,
        const std::vector<std::vector<size_t> > & snpset
) {
    int n = (int) _binary_sums(col);
    double p = _binary_probs(col);
    double score = 0.0;
    for (auto geneset = snpset.begin(); geneset != snpset.end(); geneset++) {
        int k = 0;
        for (auto id = geneset->begin(); id != geneset->end(); id++) {
            if (_expression(*id, col) > 0) {
                k++;
            }
        }
        score += -log10(gsl_ran_binomial_pdf(k, p, n));
    }
    return std::isfinite(score) ? score : 0.0;
}
