/* [X]  Read the user's SNPs.
 *      Store the names in a std::set<std::string>.
 *
 * [X]  Read the null SNPs.
 *      Store the names in a std::set<std::string>.
 *
 * [X]  Read the SNP intervals.
 *      Store the names in a std::set<std::string>.
 *      Store the intervals in a std::map<std::string, std::pair<int> >.
 *
 * [X]  Read the expression.
 *      Store the data in a Matrix.
 *      Store the names in a std::map<std::string, int> with the
 *      index of each gene name.
 *
 * [X]  Read the gene intervals.
 *      Store the intervals in a std::map<std::string, std::pair<int> >.
 *      Skip the genes that are not found in the expression map.
 * 
 * [X]  Report the user's SNP intervals and overlapping genes.
 * [X]  Put the user's absent SNPs on chromosome NA.
 *
 * [X]  Drop all SNP intervals except those in the null SNP set.
 *
 * [X]  Read the condition columns.
 *      Store the names in a std::set<std::string>.
 *
 * [X]  Report the condition columns that are absent.
 *
 * Condition the expression matrix.
 *
 * [X]  Normalize the expression matrix.
 *
 * [X]  Percentile rank the expression matrix.
 *
 * [X]  Create a std::map<int, std::vector<int> >
 *      These are bins for storing gene sets with a given number of genes.
 *
 * [X]  Create a function that generates a random SNP set (set of gene sets).
 *
 * [ ]  Create a function to calculate p-values for a single gene set.
 *
 * [ ]  Create a function to calculate p-values for a SNP set.
 *
 * [ ]  Create a function to calculate scores for each column in the expression.
 *
 * [ ]  Compare, in parallel, random SNP set scores to the user's scores.
 *
 */

#ifndef _SNPSPEC_H
#define _SNPSPEC_H

#include <Eigen/Dense>
#include "IntervalTree.h"
#include "common.h"

using namespace Eigen;

class snpspec
{
public:
    snpspec(
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
    );

    void read_names(
        std::string filename,
        std::set<std::string> & names
    );

    void read_bed_intervals(
        std::string filename,
        std::map<std::string, genomic_interval> & intervals
    );
    
    void read_gct(
            std::string filename,
            std::vector<std::string> & row_names,
            std::vector<std::string> & col_names,
            MatrixXd & data
    );

    void read_bed_interval_tree(
            std::string filename,
            const std::vector<std::string> & whitelist,
            std::map<std::string, IntervalTree<size_t> > & tree
    );

    void report_user_snp_genes(const std::string & filename);

    void drop_snp_intervals();

    void report_missing_conditions();

    void bin_genesets();

    std::vector<std::vector<size_t> > generate_snpset();

    MatrixXd geneset_pvalues_binary(std::vector<size_t> & geneset);

    double score_binary(
            const size_t & col,
            const std::vector<std::vector<size_t> > & snpset
    );

private:
    std::set<std::string>
        _user_snp_names,
        _null_snp_names,
        _snp_names,
        _condition_names;

    std::map<std::string, genomic_interval>
        _snp_intervals;

    std::map<std::string, IntervalTree<size_t> >
        _gene_interval_tree;

    MatrixXd
        _expression;

    // For binary expression, pre-calculate the column sums and those sums
    // divided by the number of rows.
    // For any expression, store the column scores for the user's SNPs.
    VectorXd
        _binary_sums,
        _binary_probs;

    // The row and column names of the provided GCT expression file.
    std::vector<std::string>
        _row_names,
        _col_names;

    // The genesets for the user's SNPs.
    std::vector<std::vector<size_t> >
        _user_genesets;

    // Each of the user's SNPs corresponds to a geneset. These are the sizes
    // of the genesets.
    std::vector<size_t>
        _user_snp_geneset_sizes;

    // Put genesets into bins, where the key to a bin is the size of the
    // contained genesets in that bin.
    std::map<int, std::vector<std::vector<size_t> > >
        _geneset_bins;

    // Is the first column of the expression filled with 1s and 0s?
    bool
        _binary_expression;
};

#endif
