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

        void report_user_snp_genes(const std::string & filename, int slop);

        void drop_snp_intervals();

        void report_missing_conditions();

        void bin_genesets(int slop, int max_genes);

        std::vector<std::vector<size_t> > generate_snpset();

        MatrixXd geneset_pvalues_binary(std::vector<size_t> & geneset);

        double score_binary(
            const size_t & col,
            const std::vector<std::vector<size_t> > & snpset
        );

        double score_quantitative(
            const size_t & col,
            const std::vector<std::vector<size_t> > & snpset
        );

        void calculate_pvalues(
            std::string filename,
            std::vector<std::vector<size_t> > genesets,
            long min_observations,
            long max_iterations,
            long replicates
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

        // Each of the user's SNPs corresponds to a geneset. These are the
        // sizes of the genesets.
        std::vector<size_t>
        _user_geneset_sizes;

        // Put genesets into bins, where the key to a bin is the size of the
        // contained genesets in that bin.
        std::map<int, std::vector<std::vector<size_t> > >
        _geneset_bins;

        // Is the first column of the expression filled with 1s and 0s?
        bool
        _binary_expression;
};

#endif
