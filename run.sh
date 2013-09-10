#!/bin/bash

time ./snpspec \
    --snps ~/dwork/snpspec/data/Hu2011/RA_SNPs.txt \
    --expression ~/dwork/snpspec/data/GO/allHumanGO.gct.gz \
    --gene-intervals ~/dwork/snpspec/reference/genes.hg19.ucsc.bed \
    --null-snps ~/dwork/snpspec/reference/TGP.pruned.txt.gz \
    --out out \
    --snp-intervals ~/dwork/snpspec/reference/TGP.bed.gz \
    --slop 250e3 \
    --condition ~/dwork/snpspec/data/GO/first_1000.txt \
    --processes 2 \
    --permutations 1e6 \
    | tee output.txt
