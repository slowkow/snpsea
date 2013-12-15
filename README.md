SNPsea
=======

SNP set enrichment analysis for condition-specificity of gene measurements or
binary annotations within trait-associated loci.

See an **[example]** of [37 SNP loci][Teslovich2010] associated with LDL
cholesterol that are significantly enriched (P=1e-4, Bonferroni=0.008) for genes
with specific expression in Liver in a test of [79 human
tissues][GSE1133].

[example]: doc/figures/Teslovich2010_Novartis2011_pvalues_barplot_25.png
[Teslovich2010]: http://www.ncbi.nlm.nih.gov/pubmed/20686565
[GSE1133]: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1133

**[Download]** the binary for 64 bit Linux.

[Download]: https://github.com/slowkow/snpsea/releases


Contents
--------

- <a href="#overview">Overview</a>
- <a href="#documentation">Documentation</a>
- <a href="#quick-start">Quick Start</a>
- <a href="#references">References</a>


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

> If trait-associated alleles impact a small number of pathogenic tissues or
> cell types in specific conditions, then the subset of genes with critical
> functions in those pathogenic cell types are likely to be within
> trait-associated genomic loci.

We assume that a gene's specificity to a given cell type or condition is
a reasonable indicator of the gene's importance to the function of that
particular cell type or condition.


Documentation
-------------

Read the manual: [HTML] or [PDF]

See <http://www.broadinstitute.org/mpg/snpsea> for more information.

[HTML]: http://www.broadinstitute.org/mpg/snpsea/SNPsea_manual.html
[PDF]: https://github.com/slowkow/snpsea/blob/master/doc/SNPsea_manual.pdf?raw=true


Quick Start
-----------

    git clone https://github.com/slowkow/snpsea.git
    cd snpsea/src
    make
    ../bin/snpsea


References
----------

Please cite the following publication if you use this method:

> Slowikowski, K. et al. *SNPsea: test trait-associated loci for enrichment of
> condition-specificity of gene measurements or binary annotations.*
> Manuscript in progress.

See additional information and examples here:

> Hu, X. et al. *Integrating autoimmune risk loci with gene-expression data
> identifies specific pathogenic immune cell subsets.* The American Journal
> of Human Genetics 89, 496–506 (2011). [PubMed][Hu2011]

[Hu2011]: http://www.ncbi.nlm.nih.gov/pubmed/21963258
