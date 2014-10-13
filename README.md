SNPsea: an algorithm to identify cell types, tissues, and pathways affected by risk loci
========================================================================================

**Home Page:** <http://www.broadinstitute.org/mpg/snpsea>

**Documentation:** [HTML] | [PDF] | [Epub]

**Executable:** [snpsea-v1.0.3.tar.gz][exec]

**Data:** [SNPsea_data_20140520.zip][data]


Citation
--------

If you benefit from this method, please cite:

> Slowikowski, K. et al. **SNPsea: an algorithm to identify cell types,
> tissues, and pathways affected by risk loci.** Bioinformatics (2014).
> doi:[10.1093/bioinformatics/btu326][Slowikowski2014]

See the first description of the algorithm and additional examples here:

> Hu, X. et al. *Integrating autoimmune risk loci with gene-expression data
> identifies specific pathogenic immune cell subsets.* The American Journal
> of Human Genetics 89, 496–506 (2011). [PubMed][Hu2011]

[Hu2011]: http://www.ncbi.nlm.nih.gov/pubmed/21963258
[Slowikowski2014]: http://bioinformatics.oxfordjournals.org/content/early/2014/05/10/bioinformatics.btu326


Description
-----------

SNPsea is an algorithm to identify cell types and pathways likely to be
affected by risk loci. It requires a list of SNP identifiers and a matrix of
genes and conditions.

Genome-wide association studies (GWAS) have discovered multiple genomic loci
associated with risk for different types of disease. SNPsea provides a simple
way to determine the types of cells influenced by genes in these risk loci.

Suppose disease-associated alleles influence a small number of pathogenic cell
types. We hypothesize that genes with critical functions in those cell types
are likely to be within risk loci for that disease. We assume that a gene's
specificity to a cell type is a reasonable indicator of its importance to the
unique function of that cell type.

First, we identify the genes in linkage disequilibrium (LD) with the given
trait-associated SNPs and score the gene set for specificity to each cell
type. Next, we define a null distribution of scores for each cell type by
sampling random SNP sets matched on the number of linked genes. Finally, we
evaluate the significance of the original gene set's specificity by comparison
to the null distributions: we calculate an exact permutation p-value.

SNPsea is a general algorithm. You may provide your own:

1. Continuous gene matrix with gene expression profiles (or other values).
2. Binary gene annotation matrix with presence/absence 1/0 values.

We provide you with three expression matrices and one annotation matrix. See
the [Data][manualdata] section of the [Manual][HTML].

The columns of the matrix may be tissues, cell types, GO annotation codes, or
other *conditions*. Continuous matrices *must* be normalized before running
SNPsea: columns must be directly comparable to each other.


License
-------

[GNU GPLv3][license]

[license]: https://github.com/slowkow/snpsea/blob/master/LICENSE

[exec]: https://github.com/slowkow/snpsea/archive/v1.0.3.tar.gz
[data]: http://files.figshare.com/1382662/SNPsea_data_20140520.zip

[HTML]: http://snpsea.readthedocs.org/en/latest/
[manualdata]: http://snpsea.readthedocs.org/en/latest/data.html
[PDF]: https://readthedocs.org/projects/snpsea/downloads/pdf/latest/
[Epub]: https://readthedocs.org/projects/snpsea/downloads/epub/latest/

