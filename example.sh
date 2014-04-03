#!/usr/bin/env bash
# example.sh
# Kamil Slowikowski
# April 3, 2014
#
# An example to demonstrate the usage of SNPsea.

options=(
    --snps              Red_blood_cell_count-Harst2012-45_SNPs.gwas
    --gene-matrix       GeneAtlas2004.gct.gz
    --gene-intervals    NCBIgenes2013.bed.gz
    --snp-intervals     TGP2011.bed.gz
    --null-snps         Lango2010.txt.gz
    --out               out
    --slop              10e3
    --threads           2
    --null-snpsets      0
    --min-observations  100
    --max-iterations    1e7
)

# Run SNPsea.
./bin/snpsea ${options[*]}

# Create a horizontal bar plot of condition p-values.
./bin/snpsea-barplot out
