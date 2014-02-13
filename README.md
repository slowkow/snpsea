SNPsea: an algorithm to identify cell types, tissues, and pathways affected by risk loci
========================================================================================

| <http://www.broadinstitute.org/mpg/snpsea> | Reference Manual [HTML] or [PDF] |
|:---:|:---:|

[HTML]: http://www.broadinstitute.org/mpg/snpsea/SNPsea_manual.html
[PDF]: https://github.com/slowkow/snpsea/blob/master/doc/SNPsea_manual.pdf?raw=true

- <a href="#citation">Citation</a>
- <a href="#quick-start">Quick Start</a>

Overview
--------

SNPsea is a general algorithm that may be used to identify cell types,
tissues, and pathways likely to be affected by risk loci.

For example, with a gene expression matrix containing expression profiles for
multiple cell types, we identify genes in linkage disequilibrium with
trait-associated SNPs and score them for specificity to each cell type.
We compare the score to a null distribution by sampling random SNP sets
matched on the number of linked genes. To evaluate significance, we calculate
an exact permutation p-value.

This implementation is generalized, so you may provide:

-   a continuous gene matrix with gene expression (or any other values)
-   a binary gene matrix with presence absence (1, 0) values.

The columns of the matrix could be tissues, cell types, GO annotation codes,
or any other types of *conditions*. Continuous matrices *must* be normalized
before running SNPsea so that columns are directly comparable to each other.

In general, this analysis is appropriate when you are interested in testing
for enrichment of condition-specificity of genes linked to a set of
trait-associated SNPs.

If trait-associated alleles impact a small number of pathogenic tissues or
cell types, we hypothesize that the subset of genes with critical functions in
those pathogenic cell types are likely to be within trait-associated loci.

We assume that a gene's specificity to a given cell type or condition is
a reasonable indicator of the gene's importance to its function.


Citation
--------

If you benefit from this method, please cite:

> Slowikowski, K. et al. *SNPsea: test trait-associated loci for enrichment of
> condition-specificity of gene measurements or binary annotations.*
> Manuscript in progress.

See additional information and examples here:

> Hu, X. et al. *Integrating autoimmune risk loci with gene-expression data
> identifies specific pathogenic immune cell subsets.* The American Journal
> of Human Genetics 89, 496–506 (2011). [PubMed][Hu2011]

[Hu2011]: http://www.ncbi.nlm.nih.gov/pubmed/21963258


Quick Start
-----------

-   **Binary for Linux**: <https://github.com/slowkow/snpsea/releases>

-   **Data**: <http://dx.doi.org/10.6084/m9.figshare.871430>

- - -

To build the source code:

1.  Install the [dependencies]:

```bash
#       Install Python libraries.
pip install docopt numpy pandas matplotlib
#       Change the graphical backend for matplotlib.
perl -i -pe 's/^(\s*(backend).*)$/#$1\n$2:Agg/' ~/.matplotlib/matplotlibrc
#       Install R libraries.
R -e 'install.packages(c("data.table", "reshape2", "gap", "ggplot2"))'
#       Install C++ libraries.
sudo apt-get install build-essential libopenmpi-dev libgsl0-dev
    #       Mac: 
    # sudo port selfupdate; sudo port install gcc48 openmpi gsl
    #       Broad Institute
    # use .gcc-4.8.1 .openmpi-1.4 .gsl-1.14
```

2.  Download and compile SNPsea:

```bash
git clone https://github.com/slowkow/snpsea.git
cd snpsea/src
make
# Run!
../bin/snpsea
```

[dependencies]: http://www.broadinstitute.org/mpg/snpsea/SNPsea_manual.html#c-libraries


