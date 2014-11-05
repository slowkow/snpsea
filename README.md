SNPsea: an algorithm to identify cell types, tissues, and pathways affected by risk loci
========================================================================================

**Home Page:** <http://www.broadinstitute.org/mpg/snpsea>

**Documentation:** [HTML] | [PDF] | [Epub]

**Executable:** [snpsea-v1.0.3.tar.gz][exec]

**Data:** [SNPsea_data_20140520.zip][data]

**License:** [GNU GPLv3][license]


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

Example
-------

![SNPsea results for RBC count-associated SNPs in the Gene Atlas.][example]

[example]: https://github.com/slowkow/snpsea/blob/master/docs/figures/Red_blood_cell_count-Harst2012-45_SNPs-GeneAtlas2004-single-pvalues_barplot.png

The heatmap shows Pearson correlation coefficients between pairs of tissue
expression profiles. The blue bars show p-values. Statistically significant
p-values cross the Bonferroni multiple testing threshold (black line).

We identified *BM-CD71+Early Erythroid* as the cell type with most significant
enrichment (P < 2e-7) for cell type-specific gene expression relative to 78
other tissues in the Gene Atlas ([Su *et al.* 2004][Su2004]).

SNPsea tested the genes in linkage disequilibrium (LD) with 45 input SNPs
associated with count of red blood cells (P <= 5e-8 in Europeans) ([Harst
*et al.* 2012][Harst2012]). For each of the 79 cell types in the Gene Atlas,
we tested a maximum of 1e7 null SNP sets where each null SNP was matched to
an input SNP on the number of genes in LD.

[Harst2012]: http://www.ncbi.nlm.nih.gov/pubmed/23222517
[Su2004]: http://www.ncbi.nlm.nih.gov/pubmed/15075390

We ran SNPsea like this:

```bash
options=(
    --snps              Red_blood_cell_count-Harst2012-45_SNPs.gwas
    --gene-matrix       GeneAtlas2004.gct.gz
    --gene-intervals    NCBIgenes2013.bed.gz
    --snp-intervals     TGP2011.bed.gz
    --null-snps         Lango2010.txt.gz
    --out               out
    --slop              10e3
    --threads           8
    --null-snpsets      0
    --min-observations  100
    --max-iterations    1e7
)
snpsea ${options[*]}

# Time elapsed: 2 minutes 36 seconds

# Create the figure shown above:
snpsea-barplot out
```

[license]: https://github.com/slowkow/snpsea/blob/master/LICENSE

[exec]: https://github.com/slowkow/snpsea/archive/v1.0.3.tar.gz
[data]: http://dx.doi.org/10.6084/m9.figshare.871430

[HTML]: http://snpsea.readthedocs.org/en/latest/
[manualdata]: http://snpsea.readthedocs.org/en/latest/data.html
[PDF]: https://readthedocs.org/projects/snpsea/downloads/pdf/latest/
[Epub]: https://readthedocs.org/projects/snpsea/downloads/epub/latest/

