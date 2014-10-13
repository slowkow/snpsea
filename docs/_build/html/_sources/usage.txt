Usage
-----

Here is a `Bash <http://www.gnu.org/software/bash/manual/bashref.html>`__
script with a usage example:

.. code:: bash

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

SNPsea will test SNPs associated with Red blood cell count for
tissue-specific expression of linked genes across 79 human tissues in
the Gene Atlas expression matrix. Each tissue will be tested up to 10
million times with matched random SNP sets, or testing will stop for a
tissue if 100 matched SNP sets achieve a higher specificity score than
the user's SNPs.

Options
~~~~~~~

All input files may optionally be compressed with
```gzip`` <http://www.gzip.org/>`__.

Required
^^^^^^^^

::


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

Optional
^^^^^^^^

::

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

Input File Formats
~~~~~~~~~~~~~~~~~~

``--snps ARG``
^^^^^^^^^^^^^^

You must provide one or more comma-separated text files. SNP identifiers
must be listed one per line. Lines starting with ``#`` are skipped. If
the file has no header, the first column is assumed to contain SNP
identifiers. Otherwise, SNPsea looks for a column named (case-sensitive)
``SNP`` or ``snp`` or ``name`` or ``marker``.

::

    head Red_blood_cell_count-Harst2012-45_SNPs.gwas

    # Harst et al. 2012
    # doi:10.1038/nature11677
    # PMID: 23222517
    # 45 SNPs associated with red blood cell count (RBC) taken from Table 1.
    # Positions are on hg19. SNPs are included if $P \le 5e-8$.
    CHR POS SNP P
    chr1    40069939    rs3916164   3e-10
    chr1    158575729   rs857684    4e-16
    chr1    199007208   rs7529925   8e-09
    chr1    248039451   rs3811444   5e-10

Instead of providing a file with SNPs, you may use "randomN" like this:

::

    --snps random20

to sample 20 random SNPs from the **``--snp-intervals``** file.

``--gene-matrix ARG``
^^^^^^^^^^^^^^^^^^^^^

You must provide a single gene matrix that must be in
`GCT <http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/gct>`__
format.

::

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

``--condition ARG`` (Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may provide column names present in the **``--gene-matrix``** file,
one per line. The matrix will be conditioned on these columns before the
analysis is performed to help you identify secondary signals independent
of these columns. Binary (0, 1) matrices will not be conditioned.

::

    head conditions.txt

    Whole_Blood

``--gene-intervals ARG``
^^^^^^^^^^^^^^^^^^^^^^^^

You must provide gene intervals in BED format with a fourth column that
contains the same gene identifiers as those present in the Name column
of the **``--gene-matrix``**
`GCT <http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/gct>`__
file. Only the first four columns are used.

::

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

``--snp-intervals ARG``
^^^^^^^^^^^^^^^^^^^^^^^

SNP linkage intervals must be specified in BED format and include a
fourth column with the SNP identifiers. The linkage intervals assigned
to the trait-associated SNPs you provide with **``--snps``** are taken
from this file.

::

    zcat TGP2011.bed.gz | head

    chr1    0   254996  rs113759966
    chr1    0   254996  rs114420996
    chr1    0   254996  rs114608975
    chr1    0   254996  rs115209712
    chr1    0   254996  rs116400033
    chr1    0   254996  rs116504101
    chr1    0   254996  rs12184306
    chr1    0   254996  rs12184307
    chr1    0   254996  rs138808727
    chr1    0   254996  rs139113303

``--null-snps ARG``
^^^^^^^^^^^^^^^^^^^

The null SNPs file must have one SNP identifier per line. Only the first
column is used. The identifiers must be a subset of the identifiers in
**``--snp-intervals``**.

::

    zcat Lango2010.txt.gz | head

    rs58108140  chr1    10583
    rs180734498 chr1    13302
    rs140337953 chr1    30923
    rs141149254 chr1    54490
    rs2462492   chr1    54676
    rs10399749  chr1    55299
    rs189727433 chr1    57952
    rs149755937 chr1    59040
    rs77573425  chr1    61989
    rs116440577 chr1    63671

Output Files
~~~~~~~~~~~~

The usage example shown above produces the following output files:

::

    out/
        args.txt
        condition_pvalues.txt
        null_pvalues.txt
        snp_condition_scores.txt
        snp_genes.txt

``args.txt``
^^^^^^^^^^^^

The command line arguments needed to reproduce the analysis.

::

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

Repeat the analysis:

::

    snpsea --args args.txt

``condition_pvalues.txt``
^^^^^^^^^^^^^^^^^^^^^^^^^

The p-values representing enrichment of condition-specificity for the
given SNPs.

::

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

``null_pvalues.txt``
^^^^^^^^^^^^^^^^^^^^

If the argument for **``--snps``** is the name of a file, the p-values
for null matched SNP sets. You can compare these null results to the
results for your trait-associated SNPs.

If the argument for **``--snps``** is "randomN" where N is some integer,
like "random20" the p-values for random unmatched SNP sets, each with N
SNPs.

The fifth column is the replicate index. The number of replicates
performed is specified with **``--null-snpsets INT``**.

::

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

``snp_genes.txt``
^^^^^^^^^^^^^^^^^

Each SNP's linkage interval and overlapping genes. If a SNP is not found
in the reference file specified with **``--snp-intervals``**, then the
name of the SNP will be listed and the other columns will contain
``NA``.

::

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

``snp_condition_scores.txt``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each SNP, condition, gene with greatest specificity to that condition,
and score for the SNP-condition pair, adjusted for the number of genes
overlapping the given SNP's linkage interval.

::

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

