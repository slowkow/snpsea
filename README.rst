=======
SNPsea
=======

SNP set enrichment analysis for condition-specificity of gene measurements or
binary annotations within trait-associated loci.

.. contents::


Overview
--------

SNPsea is a C++ program that executes a nonparametric permutation analysis for
testing enrichment of condition-specific expression of genes near
trait-associated SNPs.

The implementation described here is generalized, so you may provide
a continuous gene matrix with gene expression (or any other measurements) or
a binary gene matrix with presence absence (1, 0) values. The columns of the
matrix might be tissues, cell types, GO annotation codes, or any other types
of conditions. In general, this analysis is appropriate when you are
interested in testing for enrichment of condition-specificity of genes linked
to a given set of trait-associated SNPs.

Genome-wide association analyses have identified disease and trait loci across
the genome, thereby implicating many nearby linked genes as potentially
relevant to the trait.

The following hypothesis is tested by this analysis:

    If trait-associated alleles impact a small number of pathogenic tissues or
    cell types in specific conditions, then the subset of genes with critical
    functions in those pathogenic cell types are likely to be within
    trait-associated genomic loci.

We assume that a gene's specificity to a given cell type or condition is
a reasonable indicator of the gene's importance to the function of that
particular cell type or condition.


References
----------

Please cite the following publication if you use this method:

    Slowikowski, K et al. SNPsea: SNP set enrichment analysis for
    condition-specificity of gene measurements or binary annotations within
    trait-associated loci. Manuscript in preparation.


Quick Start
-----------

    cd src
    make
    ../bin/snpsea


Documentation
-------------

Read the manual_ and see http://www.broadinstitute.org/mpg/snpsea for more
information.

.. _manual: https://github.com/slowkow/snpsea/blob/master/doc/SNPsea_manual.pdf?raw=true
