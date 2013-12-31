// Copyright (c) 2013 Kamil Slowikowski
// See LICENSE for GPLv3 license.

#ifndef _SNPSEA_H
#define _SNPSEA_H

#define SNPSEA_VERSION "v1.0.2"

#include <Eigen/Dense>
#include "IntervalTree.h"
#include "common.h"

using namespace Eigen;

class snpsea
{
    public:
        snpsea(
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
        );

        void write_args(
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
            std::unordered_map<std::string, genomic_interval> & intervals
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
            std::unordered_map<std::string, IntervalTree<ulong> > & tree
        );

        void overlap_genes(
            std::set<std::string> & snp_names,
            std::set<std::string> & absent_snp_names,
            std::unordered_map<std::string, std::vector<ulong> > & genesets,
            std::vector<ulong> & geneset_sizes,
            ulong slop
        );

        void merge_user_snps(
            std::set<std::string> & snp_names,
            std::unordered_map<std::string, std::vector<ulong> > & genesets,
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
            const std::unordered_map<std::string, std::vector<ulong> > genesets
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
        // The set of SNPs provided by the user.
        _user_snp_names,
        // Separate the SNPs absent from --snp-intervals.
        _user_absent_snp_names,
        // Separate the SNPs whose intervals hit 0 genes.
        _user_naked_snp_names,
        // Names of SNPs from --null-snps.
        _null_snp_names,
        // Names of SNPs from --snp-intervals.
        _snp_names,
        // Names of conditions in --conditions.
        _condition_names;

        // Name of a SNP => genomic interval.
        std::unordered_map<std::string, genomic_interval>
        _snp_intervals;

        // Name of a chromosome => interval tree.
        std::unordered_map<std::string, IntervalTree<ulong> >
        _gene_interval_tree;

        MatrixXd
        _gene_matrix;

        // For a binary gene matrix, pre-calculate the column sums and those
        // sums divided by the number of rows.
        VectorXd _binary_sums;

        // For a continuous gene matrix, store the column scores for the
        // user's SNPs.
        VectorXd _binary_probs;

        // The row and column names of the provided GCT gene matrix file.
        std::vector<std::string>
        _row_names,
        _col_names;

        // The genesets for the user's SNPs.
        std::unordered_map<std::string, std::vector<ulong> >
        _user_genesets;

        // Each of the user's SNPs corresponds to a geneset. These are the
        // sizes of the genesets.
        std::vector<ulong>
        _user_geneset_sizes;

        // Put genesets into bins, where the key to a bin is the size of the
        // contained genesets in that bin.
        std::unordered_map<ulong, std::vector<std::vector<ulong> > >
        _geneset_bins;

        // Is the first column of the gene matrix filled with 1s and 0s?
        bool
        _binary_gene_matrix;
};

#endif
