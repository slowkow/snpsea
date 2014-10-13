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

<table style="margin:auto;text-align:center" cellpadding="20">
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

<table style="margin:auto;text-align:center" cellpadding="20">
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

<table style="margin:auto" cellpadding="20">
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

Supplementary Figure 1: Determining SNP linkage intervals
---------------------------------------------------------

![](figures/snp_linkage_intervals.png)

We calculated $r^{2}$ values for all pairs of SNPs within a 1 Mb sliding
window along each chromosome. Next, we assigned each of the SNPs from
The 1000 Genomes Project Phase I @consortium_integrated_2012 to a
linkage interval by identifying each SNP’s furthest upstream and
downstream neighbors with $r^{2}\ge0.5$. Finally, we extended each
interval to recombination hotspots reported by HapMap
@myers_fine-scale_2005 with recombination rate \>3 cM/Mb.

Supplementary Figure 2: Counting genes in GWAS SNP linkage intervals
--------------------------------------------------------------------

![](figures/gwas_ngenes.png)

A cumulative density plot of the number of genes overlapped by the
linkage intervals of GWAS SNPs. We downloaded the GWAS Catalog SNPs on
January 17, 2014 and selected the 11,561 SNPs present in the 1000
Genomes Project @consortium_integrated_2012. Of these SNPs, 2,119 (18%)
of them have linkage disequilibrium (LD) intervals that overlap no
genes, and 3,756 (32%) overlap a single gene. The remaining 50% of SNPs
overlap 2 or more genes. This illustrates the critical issue that many
SNPs implicate more than one gene.

Supplementary Figure 3: Choosing the $r^{2}$ threshold for linkage intervals
----------------------------------------------------------------------------

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

Supplementary Figure 4: Each trait-associated locus harbors a single associated gene
------------------------------------------------------------------------------------

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

Supplementary Figure 5: Type 1 error estimates 
-----------------------------------------------

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

<table style="margin:auto" cellpadding="20">
<tr><td>Phenotype</td><td> - - - </td><td>SNPs</td><td> - - - </td><td>Loci</td><td> - - - </td><td>Reference</td></tr>
<tr><td>Multiple sclerosis</td><td> - - - </td><td>51</td><td> - - - </td><td>47</td><td> - - - </td><td>Supp. Table A, IMSGC WTCCC 2011</td></tr>
<tr><td>Celiac disease</td><td> - - - </td><td>35</td><td> - - - </td><td>34</td><td> - - - </td><td>Table 2, Trynka, *et al.* 2011@trynka_dense_2011</td></tr>
<tr><td>HDL cholesterol</td><td> - - - </td><td>46</td><td> - - - </td><td>46</td><td> - - - </td><td>Supp. Table 2, Teslovich, *et al.* 2010 @teslovich_biological_2010</td></tr>
</table>

Supplementary Figure 6: Red blood cell count GO enrichment
----------------------------------------------------------

![](figures/Red_blood_cell_count-Harst2012-45_SNPs-GO2013-single-pvalues_barplot.png)

We observed significant enrichment for *hemopoiesis* $(2\times10^{-5})$.
The top 25 terms are shown.

Supplementary Figure 7: Multiple sclerosis
------------------------------------------

![](figures/Multiple_sclerosis-IMSGC-51_SNPs-GeneAtlas2004-single-pvalues_barplot.png)

We observed significant enrichment for 6 cell types. The top 25 of 79
are shown.

![](figures/Multiple_sclerosis-IMSGC-51_SNPs-GO2013-single-pvalues_barplot.png)

We observed significant enrichment for 52 Gene Ontology terms. The top
60 terms are shown.

Supplementary Figure 8: Celiac disease
--------------------------------------

![](figures/Celiac_disease-Trynka2011-35_SNPs-GeneAtlas2004-single-pvalues_barplot.png)

We observed significant enrichment for 3 cell types. The top 25 of 79
are shown.

![](figures/Celiac_disease-Trynka2011-35_SNPs-GO2013-single-pvalues_barplot.png)

We observed significant enrichment for 28 Gene Ontology terms. The top
40 terms are shown.

Supplementary Figure 9: HDL cholesterol
---------------------------------------

![](figures/HDL_cholesterol-Teslovich2010-46_SNPs-GeneAtlas2004-single-pvalues_barplot.png)

We observed significant enrichment for liver tissue-specific gene
expression. The top 25 of 79 are shown.

![](figures/HDL_cholesterol-Teslovich2010-46_SNPs-GO2013-single-pvalues_barplot.png)

We observed significant enrichment for 13 Gene Ontology terms. The top
25 terms are shown.
