% SNPsea Reference Manual
% Kamil Slowikowski
% February 13, 2014

\pagebreak

# Introduction

SNPsea is a general algorithm to identify cell types, tissues, and pathways
likely to be affected by risk loci.

For example, with a gene expression matrix containing expression profiles for
multiple cell types, we identify genes in linkage disequilibrium with
trait-associated SNPs and score them for specificity to each cell type.
We compare the score to a null distribution by sampling random SNP sets
matched on the number of linked genes. To evaluate significance, we calculate
an exact permutation p-value.

This implementation is generalized, so you may provide (1) a continuous gene
matrix with gene expression (or any other values) or (2) a binary gene matrix
with presence/absence 1/0 values.

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

If you benefit from this method, please cite:

> Slowikowski, K. et al. **SNPsea: an algorithm to identify cell types,
> tissues, and pathways affected by risk loci.** Manuscript in progress.

See additional examples:

> Hu, X. et al. **Integrating autoimmune risk loci with gene-expression data
> identifies specific pathogenic immune cell subsets.** The American Journal
> of Human Genetics 89, 496–506 (2011). [PubMed][Hu2011]


[Hu2011]: http://www.ncbi.nlm.nih.gov/pubmed/21963258



## Contact

Please contact me with questions and comments: <slowikow@broadinstitute.org>



## Visual Summary

### Flow Chart

![flowchart](figures/flowchart.pdf)\


This flow chart shows the input data required to perform the analysis, and
a summary of the intermediate steps.



# Installation

**Linux 64-bit**

Download the binary and run it: <https://github.com/slowkow/snpsea/releases>

**Mac**

To compile C++ code with the required dependencies, you need XCode and
MacPorts: <http://guide.macports.org/#installing.xcode>

**All other platforms**

The source code is available: <https://github.com/slowkow/snpsea>

Install the dependencies:

```bash
# Ubuntu

sudo apt-get install build-essential libopenmpi-dev libgsl0-dev

# Mac

#   Install MacPorts: http://www.macports.org/
#   Then do this:
sudo port selfupdate
sudo port install gcc48 openmpi gsl

# Broad Institute

#   You can add this to your .my.bashrc or .my.cshrc
use .gcc-4.8.1 .openmpi-1.4 .gsl-1.14
```

Download and compile the code:

```bash
#   Clone with git, so you can get updates with 'git pull'
git clone https://github.com/slowkow/snpsea.git
cd snpsea

#   or if you don't have git
curl -LOk https://github.com/slowkow/snpsea/archive/master.zip
unzip master.zip
cd snpsea-master

#   Compile.
cd src
make

#   Move the executables wherever you like.
mv ../bin/snpsea* ~/bin/
```



## Data

Download the compressed archive with data required to perform this analysis
here (138M). (The direct link to the zip shown above may be out of date,
please check.)

<http://dx.doi.org/10.6084/m9.figshare.871430>

```
GWAS SNPs

Celiac_disease-Trynka2011-35_SNPs.gwas
HDL_cholesterol-Teslovich2010-46_SNPs.gwas
Multiple_sclerosis-IMSGC-51_SNPs.gwas
Red_blood_cell_count-Harst2012-45_SNPs.gwas


GeneAtlas2004.gct.gz  # Gene Atlas 2004 gene expression matrix
GO2013.gct.gz         # Gene Ontology 2013 gene annotation matrix
ImmGen2012.gct.gz     # ImmGen 2012 gene expression matrix


NCBIgenes2013.bed.gz  # NCBI gene intervals
Lango2010.txt.gz      # LD-pruned SNPs
TGP2011.bed.gz        # 1000 Genomes Project SNP linkage intervals
```


### Celiac_disease-Trynka2011-35_SNPs.gwas

35 SNPs associated with Celiac disease taken from Table 2.
Positions are on hg19. All SNPs have P <= 5e-8.

    doi:10.1038/ng.998
    PMID: 22057235

    Trynka G, Hunt KA, Bockett NA, et al. Dense genotyping identifies and
    localizes multiple common and rare variant association signals in celiac
    disease. Nat Genet. 2011;43(12):1193-201.

    http://www.ncbi.nlm.nih.gov/pubmed/22057235


### HDL_cholesterol-Tesolvich2010-46_SNPs.gwas

46 SNPs associated with HDL taken from Supplementary Table 2.
Positions are on hg19. All SNPs have P <= 5e-8.

    doi:10.1038/nature09270
    PMID: 20686565

    Teslovich TM, Musunuru K, Smith AV, et al. Biological, clinical and
    population relevance of 95 loci for blood lipids. Nature.
    2010;466(7307):707-13.

    http://www.ncbi.nlm.nih.gov/pubmed/20686565


### Multiple_sclerosis-IMSGC-51_SNPs.gwas

51 SNPs associated with Multiple Sclerosis taken from Supplementary Table A.
Positions are on hg19. All SNPs have P <= 5e-8.

    doi:10.1038/nature10251
    PMID: 21833088

    Sawcer S, Hellenthal G, Pirinen M, et al. Genetic risk and a primary role
    for cell-mediated immune mechanisms in multiple sclerosis. Nature.
    2011;476(7359):214-9.

    http://www.ncbi.nlm.nih.gov/pubmed/21833088


### Red_blood_cell_count-Harst2012-45_SNPs.gwas

45 SNPs associated with red blood cell count (RBC) taken from Table 1.
Positions are on hg19. All SNPs have P <= 5e-8.

    doi:10.1038/nature11677
    PMID: 23222517

    Van der harst P, Zhang W, Mateo leach I, et al. Seventy-five genetic loci
    influencing the human red blood cell. Nature. 2012;492(7429):369-75.

    http://www.ncbi.nlm.nih.gov/pubmed/23222517



### GeneAtlas2004.gct.gz

Gene expression data for 79 human tissues from GSE1133. Replicates for
each tissue profile were averaged. For each gene, the single probe with
the largest minimum was selected.

    Su AI et al. A gene atlas of the mouse and human protein-encoding
    transcriptomes. Proc Natl Acad Sci U S A, 2004 Apr 9;101(16):6062-7

    http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1133


### GO2013.gct.gz

A GCT formatted gene matrix with 1s and 0s indicating presence or absence of
genes in Gene Ontology annotations. 19,111 genes in 1,751 Gene Ontology
annotations.

    http://www.geneontology.org/


### ImmGen2012.gct.gz

Gene expression data for 249 blood cell types from GSE15907. Replicates for
each cell type profile were averaged. For each gene, the single probe with
the largest minimum was selected.

    Immunological Genome Project

    http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907


### NCBIgenes2013.bed.gz

All human start and stop positions taken from:

    ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz


### Lango2010.txt.gz

A list of SNPs that span the whole genome, pruned by linkage disequilibrium
(LD). SNPsea samples null SNP sets matched on the number of genes in the
user's SNP set from this list. See this paper for more information:

    Lango allen H, Estrada K, Lettre G, et al. Hundreds of variants clustered
    in genomic loci and biological pathways affect human height. Nature.
    2010;467(7317):832-8.

    http://www.ncbi.nlm.nih.gov/pubmed/20881960


### TGP2011.bed.gz

Linkage intervals for a filtered set of SNPs from the 1000 Genomes
Project Phase 1 (May 21, 2011). SNP genotypes were obtained from the
BEAGLE release v3 website and processed to create linkage intervals for each
SNP. The linkage intervals were extended to the nearest HapMap recombination
hotspot with >3 cM/Mb recombination rate.

    http://www.1000genomes.org/
    http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes.phase1_release_v3/
    http://hapmap.ncbi.nlm.nih.gov/downloads/


[Bash]: http://www.gnu.org/software/bash/manual/bashref.html
[BEAGLE]: http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes.phase1_release_v3/
[GSE1133]: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1133
[GSE15907]: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907
[HapMap]: http://hapmap.ncbi.nlm.nih.gov/downloads/
[Lango2010]: http://www.ncbi.nlm.nih.gov/pubmed/20881960
[NCBI]: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
[Teslovich2010]: http://www.ncbi.nlm.nih.gov/pubmed/20686565
[tgp]: http://www.1000genomes.org/



## C++ Libraries

To compile SNPsea, you will need a modern C++ compiler that supports
[`c++0x`][cxx0x] and the dependencies listed below.

See: [Installation]

[intervaltree]

> a minimal C++ interval tree implementation

[Eigen]

> Eigen is a C++ template library for linear algebra: matrices, vectors,
> numerical solvers, and related algorithms.

[OpenMPI]

> MPI is a standardized API typically used for parallel and/or distributed
> computing. Open MPI is an open source, freely available implementation.

[GSL - GNU Scientific Library][gsl]

> The GNU Scientific Library (GSL) is a numerical library for C and C++
> programmers.

[GCC, the GNU Compiler][gcc]

> The GNU Compiler Collection is a compiler system produced by the GNU Project
> supporting various programming languages.

I use [`c++0x`][cxx0x] features in my C++ code, so you must use a compiler
that supports them. I compiled successfully with versions 4.6.3 (the default
version for Ubuntu 12.04) and 4.8.1.


[intervaltree]: https://github.com/slowkow/intervaltree
[Eigen]: http://eigen.tuxfamily.org
[OpenMPI]: http://www.open-mpi.org
[gsl]: http://www.gnu.org/software/gsl
[gcc]: http://gcc.gnu.org
[cxx0x]: http://gcc.gnu.org/projects/cxx0x.html



## Python Packages

To plot visualizations of the results, you will need Python 2.7 and the
packages listed below.

**Instructions:** Install with [pip]:

```
pip install docopt numpy pandas matplotlib
```

**Note:** The packages available on the Ubuntu repositories may be outdated
and might fail to work. So, avoid using `apt-get` for these dependencies.

[docopt]

> Command-line interface description language.

[numpy]

> NumPy is the fundamental package for scientific computing with Python.

[pandas]

> pandas is an open source, BSD-licensed library providing high-performance,
> easy-to-use data structures and data analysis tools for the Python
> programming language.

[matplotlib]

> matplotlib is a python 2D plotting library which produces publication
> quality figures in a variety of hardcopy formats and interactive
> environments across platforms.

**Note:** On a server with no display, please edit your [matplotlibrc] file to
use the `Agg` backend:

    perl -i -pe 's/^(\s*(backend).*)$/#$1\n$2:Agg/' ~/.matplotlib/matplotlibrc

Otherwise, you may see an error message like this:

    _tkinter.TclError: no display name and no $DISPLAY environment variable


[pip]: http://www.pip-installer.org
[docopt]: http://docopt.org/
[numpy]: http://www.numpy.org
[pandas]: http://pandas.pydata.org
[matplotlib]: http://matplotlib.org
[matplotlibrc]: http://matplotlib.org/users/customizing.html



## R Packages

Some visualizations use R and ggplot2 instead of Python and matplotlib.

**Instructions:** Start a session in R and run:

```R
install.packages(c("data.table", "reshape2", "gap", "ggplot2"))
```

[data.table]

> Extension of data.frame for fast indexing, fast ordered joins, fast
> assignment, fast grouping and list columns.

[reshape2]

> Flexibly reshape data: a reboot of the reshape package. 

[gap]

> Genetic analysis package.

[ggplot2]

> An implementation of the Grammar of Graphics.


[data.table]: http://cran.r-project.org/web/packages/data.table
[reshape2]: http://cran.r-project.org/web/packages/reshape2
[gap]: http://cran.r-project.org/web/packages/gap
[ggplot2]: http://cran.r-project.org/web/packages/ggplot2



# Usage

## Example

Here is a [Bash] script with a usage example:

```bash
options=(
    --snps              Red_blood_cell_count-Harst2012-45_SNPs.gwas
    --gene-matrix       GeneAtlas2004.gct.gz
    --gene-intervals    NCBIgenes2013.bed.gz
    --snp-intervals     TGP2011.bed.gz
    --null-snps         Lango2010.txt.gz
    --out               out
    --slop              10e3
    --threads           4
    --null-snpsets      0
    --min-observations  100
    --max-iterations    1e7
)
snpsea ${options[*]}
```

SNPs will test SNPs associated with Red blood cell count for tissue-specific
expression of linked genes across 79 human tissues in the Gene Atlas expression
matrix. Each tissue will be tested up to 10 million times with matched random
SNP sets, or testing will stop for a tissue if 100 matched SNP sets achieve a
higher specificity score than the user's SNPs.



## Options

All input files may optionally be compressed with [`gzip`][gzip].

[gzip]: http://www.gzip.org/



### Required

```

--snps ARG               Text file with SNP identifiers in the first
                         column. Instead of a file name, you may use
                         'randomN' with an integer N for a random SNP list
                         of length N.

--gene-matrix ARG        Gene matrix file in GCT format. The Name column
                         must contain the same gene identifiers as in
                         --gene-intervals.

--gene-intervals ARG     BED file with gene intervals. The fourth column
                         must contain the same gene identifiers as in
                         --gene-matrix.

--snp-intervals ARG      BED file with all known SNP intervals. The fourth
                         column must contain the same SNP identifiers as
                         in --snps and --null-snps.

--null-snps ARG          Text file with names of SNPs to sample when
                         generating null matched or random SNP sets.
                         These SNPs must be a subset of --snp-intervals.

--out ARG                Create output files in this directory. It will be
                         created if it does not already exist.
```



### Optional

```
--condition ARG          Text file with a list of columns in --gene-matrix
                         to condition on before calculating p-values. Each
                         column in --gene-matrix is projected onto each
                         column listed in this file and its projection is
                         subtracted.

--slop ARG               If a SNP interval overlaps no gene intervals,
                         extend the SNP interval this many nucleotides
                         further and try again.
                         [default: 10000]

--threads ARG            Number of threads to use.
                         [default: 1]

--null-snpsets ARG       Test this many null matched SNP sets, so you can
                         compare your results to a distribution of null
                         results.
                         [default: 0]

--min-observations ARG   Stop testing a column in --gene-matrix after
                         observing this many null SNP sets with 
                         specificity scores greater or equal to those
                         obtained with the SNP set in --snps. Increase
                         this value to obtain more accurate p-values.
                         [default: 25]

--max-iterations ARG     Maximum number of null SNP sets tested for each
                         column in --gene-matrix. Increase this value to
                         resolve smaller p-values.
                         [default: 10000]
```



## Input File Formats

### `--snps ARG`

You must provide one or more comma-separated text files. SNP identifiers must
be listed one per line. Lines starting with `#` are skipped. If the file has no
header, the first column is assumed to contain SNP identifiers. Otherwise,
SNPsea looks for a column named (case-sensitive) `SNP` or `snp` or `name` or
`marker`.

```
head Red_blood_cell_count-Harst2012-45_SNPs.gwas

# Harst et al. 2012
# doi:10.1038/nature11677
# PMID: 23222517
# 45 SNPs associated with red blood cell count (RBC) taken from Table 1.
# Positions are on hg19. SNPs are included if P <= 5e-8.
CHR	POS	SNP	P
chr1	40069939	rs3916164	3e-10
chr1	158575729	rs857684	4e-16
chr1	199007208	rs7529925	8e-09
chr1	248039451	rs3811444	5e-10
```

Instead of providing a file with SNPs, you may use "randomN" like this:

```
--snps random20
```

to sample 20 random SNPs from the **`--snp-intervals`** file.


### `--gene-matrix ARG`

You must provide a single gene matrix that must be in [GCT][gct] format.

[gct]: http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/gct


```
zcat GeneAtlas2004.gct.gz | cut -f1-4 | head

#1.2
17581  79
Name   Description  Colorectal_Adenocarcinoma  Whole_Blood
1      A1BG         115.5                      209.5
2      A2M          85                         328.5
9      NAT1         499                        1578
10     NAT2         115                        114
12     SERPINA3     419.5                      387.5
13     AADAC        125                        252.5
14     AAMP         2023                       942.5
```


### `--condition ARG` (Optional)

You may provide column names present in the **`--gene-matrix`** file, one per
line. The matrix will be conditioned on these columns before the analysis is
performed to help you identify secondary signals independent of these columns.
Binary (0, 1) matrices will not be conditioned.

```
head conditions.txt

Whole_Blood
```


### `--gene-intervals ARG`

You must provide gene intervals in BED format with a fourth column that
contains the same gene identifiers as those present in the Name column of the
**`--gene-matrix`** GCT file. Only the first four columns are used.

```
zcat NCBIgenes2013.bed.gz | head

chr1  10003485   10045555   64802      NMNAT1
chr1  100111430  100160096  54873      PALMD
chr1  100163795  100164756  100129320  HMGB3P10
chr1  100174205  100232185  391059     FRRS1
chr1  10027438   10027515   100847055  MIR5697
chr1  100308165  100308317  100270894  RPL39P9
chr1  100315632  100389578  178        AGL
chr1  100433941  100435837  730081     LOC730081
chr1  100435344  100492534  23443      SLC35A3
chr1  100503669  100548932  64645      HIAT1
```


### `--snp-intervals ARG`

SNP linkage intervals must be specified in BED format and include a fourth
column with the SNP identifiers. The linkage intervals assigned to the
trait-associated SNPs you provide with **`--snps`** are taken from this file.

```
zcat TGP2011.bed.gz | head

chr1	0	254996	rs113759966
chr1	0	254996	rs114420996
chr1	0	254996	rs114608975
chr1	0	254996	rs115209712
chr1	0	254996	rs116400033
chr1	0	254996	rs116504101
chr1	0	254996	rs12184306
chr1	0	254996	rs12184307
chr1	0	254996	rs138808727
chr1	0	254996	rs139113303
```

### `--null-snps ARG`

The null SNPs file must have one SNP identifier per line. Only the first
column is used. The identifiers must be a subset of the identifiers in
**`--snp-intervals`**.

```
zcat Lango2010.txt.gz | head

rs58108140	chr1	10583
rs180734498	chr1	13302
rs140337953	chr1	30923
rs141149254	chr1	54490
rs2462492	chr1	54676
rs10399749	chr1	55299
rs189727433	chr1	57952
rs149755937	chr1	59040
rs77573425	chr1	61989
rs116440577	chr1	63671
```



## Output Files

The usage example shown above produces the following output files:

```
out/
    args.txt
    condition_pvalues.txt
    null_pvalues.txt
    snp_condition_scores.txt
    snp_genes.txt
```


### `args.txt`

The command line arguments needed to reproduce the analysis.

```
cat args.txt

# SNPsea v1.0.2
--snps             Red_blood_cell_count-Harst2012-45_SNPs.gwas
--gene-matrix      GeneAtlas2004.gct.gz
--gene-intervals   NCBIgenes2013.bed.gz
--snp-intervals    TGP2011.bed.gz
--null-snps        Lango2010.txt.gz
--out              out
--score            single
--slop             100000
--threads          8
--null-snpsets     0
--min-observations 100
--max-iterations   10000000
```

Repeat the analysis:

```
snpsea --args args.txt
```


### `condition_pvalues.txt`

The p-values representing enrichment of condition-specificity for the given
SNPs.

```
head condition_pvalues.txt | column -t

condition                  pvalue     nulls_observed  nulls_tested
Colorectal_Adenocarcinoma  0.933555   280             300
Whole_Blood                0.521595   156             300
BM-CD33+Myeloid            0.159772   111             700
PB-CD14+Monocytes          0.103264   154             1500
PB-BDCA4+Dentritic_cells   0.0606256  187             3100
PB-CD56+NK_cells           0.194009   135             700
PB-CD4+T_cells             0.428571   128             300
PB-CD8+T_cells             0.531561   159             300
PB-CD19+B_cells            0.226819   158             700
```


### `null_pvalues.txt`

If the argument for **`--snps`** is the name of a file, the p-values for null
matched SNP sets. You can compare these null results to the results for your
trait-associated SNPs.

If the argument for **`--snps`** is "randomN" where N is some integer, like
"random20" the p-values for random unmatched SNP sets, each with N SNPs.

The fifth column is the replicate index. The number of replicates performed is
specified with **`--null-snpsets INT`**.

```
head null_pvalues.txt | column -t

ColorectalAdenocarcinoma  0.056     84   1500  0
WholeBlood                0.236667  71   300   0
BM-CD33+Myeloid           0.55      55   100   0
PB-CD14+Monocytes         0.59      59   100   0
PB-BDCA4+Dentritic_Cells  0.59      59   100   0
PB-CD56+NKCells           0.71      71   100   0
PB-CD4+Tcells             0.383333  115  300   0
PB-CD8+Tcells             0.128571  90   700   0
PB-CD19+Bcells            0.168571  118  700   0
BM-CD105+Endothelial      0.386667  116  300   0
```


### `snp_genes.txt`

Each SNP's linkage interval and overlapping genes. If a SNP is not found in
the reference file specified with **`--snp-intervals`**, then the name of the
SNP will be listed and the other columns will contain `NA`.

```
head snp_genes.txt | column -t

chrom  start      end        snp         n_genes  genes
chr4   55364224   55408999   rs218238    0        NA
chr6   139827777  139844854  rs590856    0        NA
NA     NA         NA         rs99999999  NA       NA
chr6   109505894  109651220  rs1008084   2        8763,27244
chr10  71089843   71131638   rs10159477  1        3098
chr2   111807303  111856057  rs10207392  1        55289
chr16  88831494   88903796   rs10445033  4        353,2588,9780,81620
chr7   151396253  151417368  rs10480300  1        51422
chr12  4320955    4336783    rs10849023  2        894,57103
chr15  76129642   76397903   rs11072566  4        26263,92912,123591,145957
```


### `snp_condition_scores.txt`

Each SNP, condition, gene with greatest specificity to that condition, and
score for the SNP-condition pair, adjusted for the number of genes overlapping
the given SNP's linkage interval.

```
head snp_condition_scores.txt | column -t

snp        condition                  gene   score
rs9349204  Colorectal_Adenocarcinoma  10817  0.693027
rs9349204  Whole_Blood                896    0.285864
rs9349204  BM-CD33+Myeloid            896    0.236487
rs9349204  PB-CD14+Monocytes          29964  0.340561
rs9349204  PB-BDCA4+Dentritic_cells   29964  0.411727
rs9349204  PB-CD56+NK_cells           896    0.0356897
rs9349204  PB-CD4+T_cells             896    0.38182
rs9349204  PB-CD8+T_cells             896    0.332008
rs9349204  PB-CD19+B_cells            29964  0.255196
```



## Output Visualizations


### View enrichment of tissue-specific gene expression

![A horizontal bar plot of negative log10 p-values for a test of 37
LDL-associated SNPs for enrichment of tissue-specific expression in profiles
of 79 human tissues.
](figures/Red_blood_cell_count-Harst2012-45_SNPs-GeneAtlas2004-single-pvalues_barplot.pdf)\


```bash
python bin/snpsea-barplot out
```


### View the most specifically expressed gene for each SNP-tissue pair

![A heatmap exposing the contributions of specifically expressed genes within
each SNP linkage interval to the specificity scores of each tissue.
](figures/Red_blood_cell_count-Harst2012-45_SNPs-GeneAtlas2004-single-snp_condition_heatmap.pdf)\


```bash
python bin/snpsea-heatmap out
```


### View the type 1 error rate estimates for each tissue

![A scatter plot of the observed proportion of p-values under various
thresholds after repeating the analysis with 10,000 random SNP sets.
](figures/type1error_GeneAtlas2004.pdf)\


```bash
Rscript bin/snpsea-type1error out
```



Supplementary Notes 
====================

Algorithm Details
-----------------

For a given set of SNPs associated to a single trait, SNPsea tests
whether implicated genes, in aggregate, are enriched for specificity to
an annotation or condition in a user-provided matrix of genes and
annotations/conditions. The algorithm consists of three steps:

-   **Step 1: Assigning genes to each SNP**

    -   We use linkage disequilibrium (LD) to identify the genes
        implicated by each SNP.

-   **Step 2: Calculating specificity scores**

    -   We look up implicated genes in a user-provided matrix and
        calculate a specificity score for each annotation/condition
        based on the values of these genes.

-   **Step 3: Testing significance**

    -   We compare the specificity scores to scores obtained with random
        sets of matched SNP sets and use those scores to calculate an
        empirical p-value.

### Step 1: Assigning genes to each SNP {#step-1-assigning-genes-to-each-snp .unnumbered}

Accurate analyses must address the critical issue that SNPs infrequently
implicate a locus with a single gene. Instead, SNPs implicate loci in LD
with 0 to 20 genes that may plausibly influence the associated trait
(**Supplementary Figure 2**). We determine the genes implicated by each
SNP using a previously described strategy (**Supplementary Figure 1**
and @rossin_proteins_2011). First, we define a linkage interval that
spans between the furthest neighboring SNPs within a 1 Mb window with
$r^{2}>0.5$ using European reference genomes
@consortium_integrated_2012. Next, we extend each interval to the
nearest recombination hotspots with recombination rate \>3 cM/Mb
@myers_fine-scale_2005. To address the case when no genes overlap an
interval, we provide an option to extend the interval up- and downstream
(by default 10 Kb). After each SNP has been assigned an interval and a
set of overlapping genes, we merge SNPs with shared genes into a single
locus to avoid multiple-counting of genes.

By default, we assume that one gene in each associated locus is
associated with the given trait. However, we expect some loci with
multiple genes due to poorly mapped regions in the genome, or due to
regions with high LD that span long stretches across the chromosome. We
also allow investigators the alternative to assume that all genes within
a locus are associated, and we compare results of the two alternatives
with four phenotypes (**Supplementary Figure 4**). [par:single~t~otal]

1.  The `’--score single’` method (recommended) assumes that a single
    gene in each locus is associated with the given phenotype. For each
    condition, we choose the gene in each locus with the greatest
    specificity to that condition. We assume specificity is a good way
    to represent a gene’s importance to the function of a condition*
    *(cell type or tissue).

2.  The `’--score total’` method assumes that all genes in a SNP’s
    linkage interval are associated. We account for all linked genes
    when calculating scores.

### Step 2: Calculating specificity scores {#step-2-calculating-specificity-scores .unnumbered}

SNPsea accepts matrices with continuous or binary values in each cell,
but the algorithm for a continuous matrix differs from that for a binary
matrix. Before running SNPsea, a matrix with continuous values must be
normalized so that columns are comparable. *It is not appropriate to use
this method on a “raw” matrix of expression values*, because the values
across conditions (cell types, tissues, experiments) are not directly
comparable.

#### Calculating specificity for a gene expression matrix

$$A=\text{continuous gene expression matrix}$$

We row-normalize the matrix $A$ with $m$ genes and $n$ conditions by the
L2 norm. This results in a new matrix $A'$ where each value is between 0
and 1 and represents specificity of gene $i$ in condition $j$. A value
of 1 indicates complete specificity to condition $j$ and 0 indicates
complete absence in condition $j$ as described previously
@hu_integrating_2011.

$$A'_{i,j}=\dfrac{A_{i,j}}{\sqrt{A_{i,1}^{2}+A_{i,2}^{2}+\cdots+A_{i,n}^{2}}}$$

Next, we rank each column by ascending percentile (lower rank indicates
greater specificity) to produce a new matrix $A''$. Each value is
between 0 and 1 and represents a nonparametric rank of specificity to
condition $j$ for gene $i$:

$$A''_{i,j}=\dfrac{\text{Rank}(A'_{i,j})}{m}$$

In each locus-condition pair, we choose the single gene with the highest
specificity to that condition. To correct for testing multiple genes
$g_{k}$ in each SNP locus’ set of genes $k$, we create a new matrix $P$
of probabilities. It consists of one row for each SNP locus $k$ and one
column for each condition $j$. We correct for testing multiple genes in
each SNP locus, because we choose the single gene with the greatest
specificity among $g_{k}$ genes in the locus. We also provide the option
to use all genes in each SNP locus $i\in k$.

<table style="margin:auto;margin-bottom:50px;text-align:center" cellpadding="20">
<tr>
<td>`’--score single’`</td><td>`’--score total’`</td>
</tr>
<tr>
<td>$P_{k,j}=1-(1-\text{Min}(A''_{k,j}))^{g_{k}}$</td>
<td>$P_{k,j}=-\sum_{i\in k}\text{log}A''_{i,j}$</td>
</tr>
</table>

#### Calculating specificity for a binary presence/absence matrix

$$B=\text{binary (1/0) presence/absence matrix}$$

Matrix $B$ with $m$ genes and $n$ conditions has binary values, so we do
not transform it into $B'$ or $B''$. Instead, we create a matrix $P$ of
hypergeometric probabilities. Below, $m_{j}$ genes are present in
condition $j$, $g_{k}$ genes are in locus $k$ and $g_{k,j}$ of them are
present in condition $j$. By default, we account for the presence or
absence of any of the locus’ genes in condition $j$. We also provide the
option to account for the number of genes in each locus.

<table style="margin:auto;margin-bottom:50px;text-align:center" cellpadding="20">
<tr><td></td><td>`’--score single’`</td> <td>`’--score total’`</td></tr>
<tr>
<td>$p(x)=\dfrac{{m_{j} \choose x}{m-m_{j} \choose g_{k}-x}}{{m \choose g_{k}}}$</td>
<td>$P_{k,j}=\begin{cases}
1-p(0) & g_{k,j}>0\\
1 & g_{k,j}=0
\end{cases}$</td><td>$P_{k,j}=\begin{cases}
1-\sum_{x=0}^{g_{k,j}-1}p(x) & g_{k,j}>0\\
1 & g_{k,j}=0
\end{cases}$</td>
</tr>
</table>

Finally, we define a specificity score $S_{j}$ for each condition $j$ as
the aggregate of all loci:

$$S_{j}=-\sum_{k}\text{log}P_{k,j}$$

### Step 3: Testing significance  {#step-3-testing-significance .unnumbered}

#### Analytical p-values

We previously found that the use of analytical distributions such as the
Poisson distribution can result in inaccurate p-values
@hu_integrating_2011. SNP linkage intervals, gene densities, gene sizes
and gene functions are correlated across the genome and are challenging
to model analytically. This motivated us to use a sampling strategy
instead.

#### Permutation p-values

For each condition, we use a sampling approach to calculate an empirical
p-value: the tail probability of observing a condition-specificity score
$S_{j}$ in a distribution of scores obtained with null SNP sets. We
compute specificity scores $S$ for random SNP sets where each SNP is
matched to a SNP in the user’s set on the number of linked genes. To
adequately cover genes across the genome, we sample SNP sets from a list
of LD-pruned SNPs that is a subset of the SNPs from the 1000 Genomes
Project @lango_allen_hundreds_2010. For each condition $j$, we calculate
an exact permutation p-value @phipson_permutation_2010 where $a_{j}$ is
the number of sampled SNP sets and $b_{j}$ of them have a specificity
score greater than or equal to the user’s score $S_{j}$:

$$p_{j}=\dfrac{b_{j}+1}{a_{j}+1}$$

We implemented adaptive sampling to calculate p-values more efficiently.
As each condition is tested, we perform more iterations to increase
resolution of significant p-values and fewer iterations for less
significant p-values to save computation time. Two options allow the
user to control the adaptive sampling:

1.  `’--max-iterations N’` The maximum number of iterations for each
    condition. We stop testing a condition after sampling $N$ SNP sets.

2.  `’--min-observations N’` The minimum number of observed higher
    scores $S_{j}$ required to stop sampling SNP sets for a condition
    $j$.

Data
----

Please find the data required to reproduce this analysis here:
<http://dx.doi.org/10.6084/m9.figshare.871430>

### Gene Atlas gene expression matrix {#gene-atlas-gene-expression-matrix .unnumbered}

We downloaded the data from BioGPS:
<http://plugins.biogps.org/download/gnf1h-gcrma.zip>

We averaged the expression values for tissue replicates. For each gene,
we selected the single probe with the largest minimum value. Finally, we
converted the file to [GCT
format](http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/gct).

### Gene Ontology binary presence/absence matrix {#gene-ontology-binary-presenceabsence-matrix .unnumbered}

We downloaded the OBO file from Gene Ontology (data-version: 2013-06-29,
CVS revision: 9700):

<http://www.geneontology.org>

For each gene, we climbed the hierarchy of ontology terms and applied
parental terms. If a gene is annotated with some term $T$, we also add
all of the terms that are parents of $T$. We copy terms between
homologous genes using Homologene data
(<http://www.ncbi.nlm.nih.gov/homologene>). If a mouse gene is annotated
with some term and the human homolog is not, then we copy the term to
the human gene. We discard all GO terms assigned to fewer than 100 or to
more than 1000 genes. This leaves us with a matrix of 19,111 genes and
1,751 terms.

### 1000 Genomes Project {#genomes-project .unnumbered}

We downloaded a filtered (diallelic and 5 or more copies of the minor
allele) set of markers from the BEAGLE website and calculated pairwise
LD (EUR) for all SNPs in a 1 Mb sliding window:

<http://bochet.gcc.biostat.washington.edu/beagle>

Commands
--------

<table style="margin:auto;margin-bottom:50px" cellpadding="20">
<header><td>Analysis</td><td></td><td></td><td>Figures</td></header>
<tr><td>
```bash
snpsea
--snps HDL_Teslovich2010.txt
--gene-matrix GeneAtlas2004.gct.gz
--gene-intervals NCBIgenes2013.bed.gz
--snp-intervals TGP2011.bed.gz
--null-snps Lango2010.txt.gz
--score single
--out <out>
--slop 10e3
--threads 2
--null-snpsets 0
--min-observations 100
--max-iterations 1e6
```
</td>
<td></td><td> - - - </td>
<td>
Horizontal barplot of p-values:</br></br>
`snpsea-barplot <out>`</br></br>
Heatmap of conditions and loci:</br></br>
`snpsea-heatmap <out>`</br></br>
Type 1 error rates for each condition:</br></br>
`snpsea-type1error <out>`</br>
<td></tr>
</table>

Supplementary Figures
=====================

**Supplementary Figure 1: Determining SNP linkage intervals**

![](figures/snp_linkage_intervals.png)

We calculated $r^{2}$ values for all pairs of SNPs within a 1 Mb sliding
window along each chromosome. Next, we assigned each of the SNPs from
The 1000 Genomes Project Phase I @consortium_integrated_2012 to a
linkage interval by identifying each SNP’s furthest upstream and
downstream neighbors with $r^{2}\ge0.5$. Finally, we extended each
interval to recombination hotspots reported by HapMap
@myers_fine-scale_2005 with recombination rate \>3 cM/Mb.

**Supplementary Figure 2: Counting genes in GWAS SNP linkage intervals**

![](figures/gwas_ngenes.png)

A cumulative density plot of the number of genes overlapped by the
linkage intervals of GWAS SNPs. We downloaded the GWAS Catalog SNPs on
January 17, 2014 and selected the 11,561 SNPs present in the 1000
Genomes Project @consortium_integrated_2012. Of these SNPs, 2,119 (18%)
of them have linkage disequilibrium (LD) intervals that overlap no
genes, and 3,756 (32%) overlap a single gene. The remaining 50% of SNPs
overlap 2 or more genes. This illustrates the critical issue that many
SNPs implicate more than one gene.

**Supplementary Figure 3: Choosing the $r^{2}$ threshold for linkage intervals**

We chose to use $r^{2}\geq0.5$ due to previous experience
@rossin_proteins_2011. To investigate if this choice influences SNPsea
results, we repeated the analysis of 45 red blood cell count-associated
SNPs @van_der_harst_seventy-five_2012 using 5 different thresholds
($r^{2}\ge$ 0.2, 0.4, 0.6, 0.8, 1.0). We also did this for SNPs
associated with multiple sclerosis, celiac disease, and HDL cholesterol.

<table><tr><td>
![](figures/r2_thresholds-GeneAtlas2004.png)
</td><td>
![](figures/r2_thresholds-GO2013.png)
</td></tr></table>

Gene Atlas and Gene Ontology (left and right). Each subplot has
$-\text{log}_{10}P$ for $r^{2}=1$ on the x-axis and $-\text{log}_{10}P$
on the y-axis for the $r^{2}$ threshold marked above. Grey lines are
significance thresholds after correction testing multiple conditions
(cell types, GO annotations). Black points are significant and grey are
not. We used the `’--score single’` option. Red blood cell count SNPs
are enriched for *hemopoiesis* (GO:0030097) ($P=2\times10^{-5}$) for
linkage intervals with $r^{2}=(0.6,0.8,1.0)$. This result falls below
the multiple testing threshold at $r^{2}\ge0.4$, but remains significant
at $r^{2}\ge0.5$ (see main text).

**Supplementary Figure 4: Each trait-associated locus harbors a single associated gene**

![](figures/single_vs_total-RBC-MS-CD-HDL.png)

Gene Atlas @su_gene_2004 and Gene Ontology (top and bottom). The x and y
axes are$-\text{log}_{10}P$ for `’--score single’ and ’--score total’`
SNPsea options.

We tested our assumption that a single gene in each trait-associated
locus is driving the condition-specificity signal by running SNPsea with
the `’single’` and `’total’` methods described . We lose power to detect
gene expression enrichment when we relax the single gene assumption and
instead assume that *all* genes implicated by LD are driving the signal.

This is in contrast with our observations for GO term enrichment, where
p-values appear similar between methods: `’single’` where we check for
the presence of a one or more genes in the locus anotated with a term
and `’total’` where we record the number of genes in the locus annotated
with a term.

**Supplementary Figure 5: Type 1 error estimates **

![](figures/type1error_GeneAtlas2004.png)

![](figures/type1error_GO2013.png)

We sampled 10,000 sets of 100 SNPs uniformly from a list of LD-pruned
SNPs @lango_allen_hundreds_2010. We tested each of the 10,000 sets for
enrichment of tissue-specific expression in the Gene Atlas @su_gene_2004
gene expression matrix (top) and for enrichment of annotation with Gene
Ontology terms (bottom). For each condition, we show the proportion of
the 10,000 enrichment p-values that are below the given thresholds. We
observe that the p-values are near the expected values, so the type 1
(false positive) error rate is well-calibrated.

Additional Examples
===================

We tested SNPsea with the three additional phenotypes listed below with
genome-wide significant SNPs $(P\leq5\times10^{-8})$. When multiple SNPs
implicated the same genes, we merged them into a single locus. We tested
each phenotype with the Gene Atlas and GO matrices with the
`’--score single’` option. The adjacent heatmaps show Pearson
correlation coefficients for all pairs of conditions.

<table style="margin:auto;margin-bottom:50px" cellpadding="20">
<tr><td>Phenotype</td><td> - - - </td><td>SNPs</td><td> - - - </td><td>Loci</td><td> - - - </td><td>Reference</td></tr>
<tr><td>Multiple sclerosis</td><td> - - - </td><td>51</td><td> - - - </td><td>47</td><td> - - - </td><td>Supp. Table A, IMSGC WTCCC 2011</td></tr>
<tr><td>Celiac disease</td><td> - - - </td><td>35</td><td> - - - </td><td>34</td><td> - - - </td><td>Table 2, Trynka, *et al.* 2011@trynka_dense_2011</td></tr>
<tr><td>HDL cholesterol</td><td> - - - </td><td>46</td><td> - - - </td><td>46</td><td> - - - </td><td>Supp. Table 2, Teslovich, *et al.* 2010 @teslovich_biological_2010</td></tr>
</table>

**Supplementary Figure 6: Red blood cell count GO enrichment**

![](figures/Red_blood_cell_count-Harst2012-45_SNPs-GO2013-single-pvalues_barplot.png)

We observed significant enrichment for *hemopoiesis* $(2\times10^{-5})$.
The top 25 terms are shown.

**Supplementary Figure 7: Multiple sclerosis**

![](figures/Multiple_sclerosis-IMSGC-51_SNPs-GeneAtlas2004-single-pvalues_barplot.png)

We observed significant enrichment for 6 cell types. The top 25 of 79
are shown.

![](figures/Multiple_sclerosis-IMSGC-51_SNPs-GO2013-single-pvalues_barplot.png)

We observed significant enrichment for 52 Gene Ontology terms. The top
60 terms are shown.

**Supplementary Figure 8: Celiac disease**

![](figures/Celiac_disease-Trynka2011-35_SNPs-GeneAtlas2004-single-pvalues_barplot.png)

We observed significant enrichment for 3 cell types. The top 25 of 79
are shown.

![](figures/Celiac_disease-Trynka2011-35_SNPs-GO2013-single-pvalues_barplot.png)

We observed significant enrichment for 28 Gene Ontology terms. The top
40 terms are shown.

**Supplementary Figure 9: HDL cholesterol**

![](figures/HDL_cholesterol-Teslovich2010-46_SNPs-GeneAtlas2004-single-pvalues_barplot.png)

We observed significant enrichment for liver tissue-specific gene
expression. The top 25 of 79 are shown.

![](figures/HDL_cholesterol-Teslovich2010-46_SNPs-GO2013-single-pvalues_barplot.png)

We observed significant enrichment for 13 Gene Ontology terms. The top
25 terms are shown.
