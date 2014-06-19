SNPsea: an algorithm to identify cell types, tissues, and pathways affected by risk loci
========================================================================================

| [broadinstitute.org/mpg/snpsea][broad] | [Reference Manual][manual] |
|:---:|:---:|
| [Example](#example) | [Citation](#citation) |
| [Quick Start](#quick-start) | [Description](#description) |
| Executable | Data |
| [snpsea-v1.0.3.tar.gz][exec] | [SNPsea_data_20140520.zip][data]  |

[broad]: http://www.broadinstitute.org/mpg/snpsea
[manual]: http://www.broadinstitute.org/mpg/snpsea/SNPsea_manual.html
[exec]: https://github.com/slowkow/snpsea/archive/v1.0.3.tar.gz
[data]: http://files.figshare.com/1382662/SNPsea_data_20140520.zip


Cartoon
-------

![Cartoon of SNPsea algorithm.][fig1]

[fig1]: https://raw.github.com/slowkow/snpsea/master/doc/figures/cartoon.png

This cartoon illustrates the key ideas of the algorithm:

**A|**  Step 1. Each SNP in a set of disease-associated SNPs is in linkage
        disequilibrium (LD) with multiple genes. The genes are scored, in
        aggregate, for specificity to each tissue.

**B|**  Step 2: The algorithm is repeated with random null SNP sets that are
        not associated with any phenotype. These have been selected from an
        LD-pruned list of SNPs, so the whole genome is covered.

**C|**  Step 3: The random SNP set scores form the null distributions which
        allows us to determine statistical significance for enrichment of
        specificity to a particular tissue/cell-type/condition.


Example
-------

![Example of SNPsea results.][fig2]

[fig2]: https://raw.github.com/slowkow/snpsea/master/doc/figures/Red_blood_cell_count-Harst2012-45_SNPs-GeneAtlas2004-single-pvalues_barplot.png

We identified *BM-CD71+Early Erythroid* as the cell type with most significant
enrichment (P < 2e-7) for cell type-specific gene expression relative to 78
other tissues in the Gene Atlas ([Su *et al.* 2004][Su2004]).

We used SNPsea to test the genes in linkage disequilibrium (LD) with 45 SNPs
associated with red blood cell count (P <= 5e-8) in GWAS of Europeans ([Harst
*et al.* 2012][Harst2012]). For each cell type, we tested a maximum of 1e7
null SNP sets by sampling random LD pruned SNP sets and scoring them. Each
set contains random SNPs matched to the input SNPs on the number of genes in
LD.

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


Quick Start
-----------

| Executable | Data |
|:---:|:---:|
| [snpsea-v1.0.3.tar.gz][exec] | [SNPsea_data_20140520.zip][data]  |


On Linux 64-bit, get started right away:

```bash
curl -LOk https://github.com/slowkow/snpsea/archive/v1.0.3.tar.gz
tar xf v1.0.3.tar.gz
cd snpsea-1.0.3
curl -LOk http://files.figshare.com/1382662/SNPsea_data_20140520.zip
unzip SNPsea_data_20140520.zip
bash example.sh
```

- - -

On other platforms, build the source code.

1. Install the [dependencies]:

    ```bash
    #       Install Python libraries.
    pip install docopt numpy pandas matplotlib
    #       Change the graphical backend for matplotlib.
    perl -i -pe 's/^(\s*(backend).*)$/#$1\n$2:Agg/' ~/.matplotlib/matplotlibrc
    
    #       Install R libraries.
    R -e 'install.packages(c("data.table", "reshape2", "gap", "ggplot2"))'
    
    #       Install C++ libraries.
    #       Ubuntu
    sudo apt-get update; sudo apt-get install build-essential libopenmpi-dev libgsl0-dev
    #       Mac
    sudo port selfupdate; sudo port install gcc48 openmpi gsl
    #       Broad Institute
    use .gcc-4.8.1 .openmpi-1.4 .gsl-1.14
    ```

2. Download and compile SNPsea:

    ```bash
    git clone https://github.com/slowkow/snpsea.git
    cd snpsea/src
    make
    ```

[dependencies]: http://www.broadinstitute.org/mpg/snpsea/SNPsea_manual.html#c-libraries


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
the [Data][manualdata] section of the [Reference Manual][manual].

[manualdata]: http://www.broadinstitute.org/mpg/snpsea/SNPsea_manual.html#data

The columns of the matrix may be tissues, cell types, GO annotation codes, or
other *conditions*. Continuous matrices *must* be normalized before running
SNPsea: columns must be directly comparable to each other.


License
-------

[GNU GPLv3][license]

[license]: https://github.com/slowkow/snpsea/blob/master/LICENSE
