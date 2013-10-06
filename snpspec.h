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
            ulong slop,
            int threads,
            ulong null_snpset_replicates,
            ulong min_observations,
            ulong max_iterations 
        );

        void write_args(
            std::vector<std::string> user_snpset_files,
            std::string expression_file,
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
            int which_snpset_file,
            std::ostream & stream
        );

        void read_names(
            std::string filename,
            std::set<std::string> & names
        );

        std::vector<ulong> snp_geneset(std::string, ulong slop);

        void random_snps(
            std::string filename,
            std::set<std::string> & names,
            ulong slop
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
            const std::vector<std::string> & row_names,
            std::map<std::string, IntervalTree<ulong> > & tree
        );

        void overlap_genes(
            std::set<std::string> & snp_names,
            std::set<std::string> & absent_snp_names,
            std::map<std::string, std::vector<ulong> > & genesets,
            std::vector<ulong> & geneset_sizes,
            ulong slop
        );

        void merge_user_snps(
            std::set<std::string> & snp_names,
            std::map<std::string, std::vector<ulong> > & genesets,
            std::vector<ulong> & geneset_sizes
        );

        void report_user_snp_genes(const std::string & filename);

        void drop_snp_intervals();

        void report_missing_conditions();

        void condition(
            MatrixXd & matrix,
            std::set<std::string> & col_names
        );

        void bin_genesets(ulong slop, ulong max_genes);

        std::vector<std::vector<ulong> > matched_genesets();

        std::vector<std::vector<ulong> > random_genesets(int n, ulong slop);

        MatrixXd geneset_pvalues_binary(std::vector<ulong> & geneset);

        double score_binary(
            const ulong & col,
            const std::vector<std::vector<ulong> > & genesets
        );

        double score_quantitative(
            const ulong & col,
            const std::vector<std::vector<ulong> > & genesets
        );

        void report_pvalues(
            const std::string filename,
            const std::map<std::string, std::vector<ulong> > genesets
        );

        void calculate_pvalues(
            std::string filename,
            std::vector<std::vector<ulong> > genesets,
            long min_observations,
            long max_iterations,
            long replicates
        );

    private:
        std::set<std::string>
        _user_snp_names,
        _user_absent_snp_names,
        _null_snp_names,
        _snp_names,
        _condition_names;

        // Name of a SNP => genomic interval.
        std::map<std::string, genomic_interval>
        _snp_intervals;

        // Name of a chromosome => interval tree.
        std::map<std::string, IntervalTree<ulong> >
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
        std::map<std::string, std::vector<ulong> >
        _user_genesets;

        // Each of the user's SNPs corresponds to a geneset. These are the
        // sizes of the genesets.
        std::vector<ulong>
        _user_geneset_sizes;

        // Put genesets into bins, where the key to a bin is the size of the
        // contained genesets in that bin.
        std::map<ulong, std::vector<std::vector<ulong> > >
        _geneset_bins;

        // Is the first column of the expression filled with 1s and 0s?
        bool
        _binary_expression;
};

#endif
