Algorithm Details
-----------------

SNPsea tests if genes implicated by risk loci (e.g., those discovered
through genome-wide association (GWA) studies) are specifically
expressed in some conditions over others, and if this specificity is
statistically significant. The program requires two inputs:

1. A list of SNP identifiers: rs123, 12:456, ...

2. A matrix of genes and conditions, such as:

   -  Gene expression profiles of multiple different cell types.

   -  Ontology terms and presence/absence 1/0 values for each gene in
      each term.

For example, SNPsea can be used to find tissues or cell types whose
function is likely to be influenced by genes in risk loci. If the genes
in risk loci are used in relatively few cell types, we hypothesize that
they are likely to affect those cell types’ unique functions. This
assumes that expression specificity is a good indicator of a gene’s
importance to the unique function of the cell type.

For a given set of SNPs associated to some phenotype, SNPsea tests
whether all implicated genes, in aggregate, are enriched for specificity
to a condition in a user-provided matrix of genes and
conditions/annotations. The algorithm consists of three steps:

-  **Step 1: Assigning genes to each SNP**

   -  We use linkage disequilibrium (LD) to identify the genes
      implicated by each SNP.

-  **Step 2: Calculating specificity scores**

   -  We look up implicated genes in a user-provided matrix and
      calculate a specificity score for each annotation/condition based
      on the values of these genes.

-  **Step 3: Testing significance**

   -  We compare the specificity scores to a null distribution of scores
      obtained with random sets of matched SNP sets and compute an
      empirical :math:`P`-value.

Step 1: Assigning genes to each SNP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Accurate analyses must address the critical issue that SNPs frequently
implicate a region with multiple different genes ( :ref:`Supplementary Figure
2 <fig-s2>` ). The challenge is to find evidence to show which of those genes
are associated with a given trait.

We determine the genes plausibly implicated by each trait-associated SNP using
a previously described strategy ( :ref:`Supplementary Figure 1 <fig-s1>` and
Rossin *et al.* 2011). First, we define the linkage interval for a given SNP
as the span between the furthest correlated SNPs :math:`r^{2}>0.5` (EUR)
within a 1 Mb window (1000 Genomes Consortium 2012). Next, we extend the
interval to the nearest recombination hotspots with recombination rate >3
cM/Mb (Myers *et al.* 2005). To address the case when no genes overlap an
interval, we provide an option for SNPsea to extend the interval up- and
downstream (by default 10 Kb).

Most frequently, we find multiple genes :math:`(m_{k}>1)` in a single
SNP locus :math:`k`. We expect many loci with multiple genes because of
regions with high LD across long stretches of a chromosome. Less
frequently, a locus has a single gene :math:`(m_{k}=1)`, and loci with
no genes :math:`(m_{k}=0)` are discarded.

After each SNP has been assigned an interval and a set of genes
overlapping the interval, we merge SNPs with shared genes into a single
locus to avoid multiple-counting of genes.

Two score options
^^^^^^^^^^^^^^^^^

By default, SNPsea assumes one gene in each associated locus is
associated with the given trait. We also include the option to assume
all genes within a locus are associated. We compare results of the two
options with four phenotypes ( :ref:`Supplementary Figure 4 <fig-s4>` ).

1. The ``’--score single’`` method (default option) assumes that a
   single gene in each locus is associated with the given phenotype. For
   each condition, we choose the gene in each locus with the greatest
   specificity to that condition.

2. The ``’--score total’`` method assumes that all genes in a SNP’s
   linkage interval are associated. We account for all linked genes when
   calculating scores.

Step 2: Calculating specificity scores
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SNPsea uses different algorithms for matrices with continuous or binary
values. Before running SNPsea, a matrix with continuous values must be
normalized so that columns are directly comparable. *It is not
appropriate to use this method on a “raw” matrix of expression values*.

Specificity for a matrix of continuous values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We extend an approach we have previoulsy described in detail (Hu *et
al.* 2011). Let :math:`A` denote a continuous gene expression matrix
with :math:`m` genes and :math:`n` conditions. First, we normalize the
expression of each gene by dividing each value by the L2 norm of the
genes values in different conditions.

.. math:: A'_{i,j}=\dfrac{A_{i,j}}{\sqrt{A_{i,1}^{2}+A_{i,2}^{2}+\cdots+A_{i,n}^{2}}}

The resulting matrix :math:`A'` has values :math:`A'_{i,j}` between 0
and 1 indicating specificity of gene :math:`i` to condition :math:`j`. A
value :math:`A'_{i,j}=1` indicates that gene :math:`i` is exclusively
expressed in condition :math:`j`, and :math:`A'_{i,j}=0` indicates that
gene :math:`i` is not expressed in condition :math:`j`.

Next, we transform :math:`A'` to a matrix :math:`A''` of non-parametric
condition-specificity percentiles as follows. For each condition
:math:`j`, we rank the values of :math:`A'_{,j}` in ascending order and
divide them by the number of genes :math:`m`, resulting in percentiles
between 0 and 1 where a lower value indicates greater specificity to the
given condition.

.. math:: A''_{i,j}=\dfrac{\text{Rank}_{j}(A'_{i,j})}{m}

Locus scores for a matrix of continuous values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We create a new matrix :math:`P`, where each value :math:`P_{k,j}` is a
score for a SNP locus :math:`k` and a condition :math:`j`. The locus
scores :math:`P_{,j}` for a single condition :math:`j` are approximately
uniformly distributed for a set of randomly selected loci under the
following assumption: for the set of genes in a given SNP locus
:math:`I_{k}`, the values :math:`A''_{i\in I_{k},j}` are random,
independent, and approximately uniformly distributed. We’ll come back to
this assumption later when testing significance in Step 3 below.

``’--score single’ (default)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This approach assumes one gene in each SNP locus is associated with the
trait.

For each locus-condition pair :math:`(k,j)`, we choose the single gene
:math:`i` in locus :math:`k` with greatest specificity to condition
:math:`j` among the :math:`m_{k}` genes in the locus, as previously
described in Hu *et al.* (Hu *et al.* 2011). Let :math:`g_{k}` denote
this most specific gene, so that
:math:`A''_{g_{k},j}=\text{Min}_{i\in I_{k}}(A''_{i,j})` where
:math:`I_{k}` denotes the set of genes in locus :math:`k`. If we assume
values of :math:`A''_{i\in I_{k},j}` are uniformly distributed for a
given condition :math:`j` and genes :math:`i\in I_{k}`, then the
probability of obtaining a value equal to or less than
:math:`A''_{g_{k},j}` is as follows:

.. math:: P_{k,j}=1-(1-\text{Min}_{i\in I_{k}}(A''_{i,j}))^{m_{k}}

``’--score total’``
^^^^^^^^^^^^^^^^^^^

This assumes all genes in a given SNP locus are associated with a trait — we
consider this model to be unlikely in most situations. We compute the
probability of observing values :math:`A''_{i\in I_{k}}` for some locus
:math:`k` as the product of percentiles. This assumes :math:`A''_{i\in I_{k}}`
values are uniformly distributed.

.. math::

   \begin{aligned}
   P_{k,j} & = & \int_{x}^{\infty}\Gamma(m_{k},1)\ \ \text{for}\ \ x=\sum_{i\in I_{k}}-\text{ln}A''_{i,j}\end{aligned}

Locus scores for a matrix of binary values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let :math:`B` denote a binary matrix (1=present, 0=absent) with
:math:`m` genes and :math:`n` conditions. Let :math:`m_{j}` denote the
number of genes present in condition :math:`j`. Let :math:`m_{k}` denote
the number of genes in locus :math:`k` and :math:`m_{k,j}\le m_{k}`
denote the number of genes in locus :math:`k` that are present in
condition :math:`j`.

We provide two options to calculate locus scores. By default, we account
for presence or absence of any of the :math:`m_{k}` genes in condition
:math:`j`, as shown below (``’--score single’``). Alternatively, we
account for the number of genes in a given locus (``’--score total’``).

+--------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------+
| \`’--score single’                                                             | ``’--score total’``                                                                                  |
+================================================================================+======================================================================================================+
| :math:`P_{k,j}=\begin{cases} 1-p(0) & m_{k,j}>0\\ 1 & m_{k,j}=0 \end{cases}`   | :math:`P_{k,j}=\begin{cases} 1-\sum_{x=0}^{m_{k,j}-1}p(x) & m_{k,j}>0\\ 1 & m_{k,j}=0 \end{cases}`   |
+--------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------+

where

.. math:: p(x)=\dfrac{{m_{j} \choose x}{m-m_{j} \choose m_{k}-x}}{{m \choose m_{k}}}

Condition specificity scores
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For both continuous and binary matrices, we define a specificity score
:math:`S_{j}` for each condition :math:`j` as the aggregate of
:math:`P_{k,j}` values across SNP loci:

.. math:: S_{j}=\sum_{k}-\text{log}P_{k,j}

Step 3: Testing significance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analytical p-values
^^^^^^^^^^^^^^^^^^^

We previously found that aggregating the :math:`P_{k,j}` scores and
determining a :math:`P`-value analytically from a distribution results
in inaccurate p-values (Hu *et al.* 2011). :math:`A''_{i,j}` values may
be relatively uniform genome-wide, but proximate genes often have shared
functions. The genome has a complex correlation structure of linkage
disequilibrium, gene density, gene size and function that is challenging
to model analytically. We use the sampling strategy described below
instead.

Permutation p-values
^^^^^^^^^^^^^^^^^^^^

For each condition, we use a sampling approach to calculate an empirical
p-value. This is the tail probability of observing a
condition-specificity score greater or equal to :math:`S_{j}`. We obtain
the distribution empirically with null SNP sets.

We compute specificity scores :math:`S` for random SNP sets. Each SNP in
a null set is matched to a SNP in the user’s set on the number of linked
genes. To adequately sample genes from the entire genome, we sample SNP
sets from a list of LD-pruned SNPs (subset of SNPs in 1000 Genomes
Project) (Lango Allen *et al.* 2010).

For each condition :math:`j`, we calculate an exact permutation p-value
(Phipson *et al.* 2010). Let :math:`a_{j}` denote the number of sampled
SNP sets (e.g. 10,000) and let :math:`b_{j}` denote how many null
specificity scores are greater than or equal to the user’s score
:math:`S_{j}`:

.. math:: p_{j}=\dfrac{b_{j}+1}{a_{j}+1}

We implemented adaptive sampling to calculate p-values efficiently. As
each condition is tested for significance, we increase the number of
iterations to resolve significant p-values and save computation by using
fewer iterations for less significant p-values. Two options allow the
user to control the adaptive sampling:

1. ``’--max-iterations N’`` The maximum number of iterations for each
   condition. We stop testing a condition after sampling :math:`N` SNP
   sets.

2. ``’--min-observations N’`` The minimum number of observed null
   specificity scores greater than or equal to :math:`S_{j}` required to
   stop sampling SNP sets for a condition :math:`j`.

Example
~~~~~~~

Suppose we have a gene expression matrix :math:`A`:

.. code-block:: r

   > A1 = read.table(text = "
   2.55 0.05 3.28 1.11
   2.63 4.53 4.66 3.89
   0.61 3.31 2.49 4.59
   0.82 1.27 4.47 2.31
   4.91 1.23 0.51 0.95")
   > A1
       V1   V2   V3   V4
   1 2.55 0.05 3.28 1.11
   2 2.63 4.53 4.66 3.89
   3 0.61 3.31 2.49 4.59
   4 0.82 1.27 4.47 2.31
   5 4.91 1.23 0.51 0.95

Compute the specificity (L2 norm) of each gene (row) to each condition
(column):

.. code-block:: r

   > A2 = t(apply(A1, 1, function(row) row / sqrt( sum(row ^ 2) )))
   > A2
                V1         V2         V3        V4
   [1,] 0.59293508 0.01162618 0.76267727 0.2581012
   [2,] 0.32801918 0.56499121 0.58120508 0.4851690
   [3,] 0.09818755 0.53278820 0.40079837 0.7388211
   [4,] 0.15607783 0.24173030 0.85081451 0.4396827
   [5,] 0.94873958 0.23766796 0.09854525 0.1835647

Rank the genes in each condition and convert to percentiles:

.. code-block:: r

   A3 = apply(A2, 2, function(col) rank(-col) / length(col))
   > A3
         V1  V2  V3  V4
   [1,] 0.4 1.0 0.4 0.8
   [2,] 0.6 0.2 0.6 0.4
   [3,] 1.0 0.4 0.8 0.2
   [4,] 0.8 0.6 0.2 0.6
   [5,] 0.2 0.8 1.0 1.0

Notice that gene 3 has the greatest specificity (0.74) to condition V4, so it
is assigned the lowest percentile rank (0.2).

Compute the locus scores for a SNP locus :math:`k` that overlaps genes 2 and
4, assuming that a single gene (either 2 or 4 but not both) is associated with
the trait:

.. code-block:: r

   > genes = c(2, 4)
   > P = apply(A3[genes, ], 2, function(col) 1 - (1 - min(col)) ^ length(col))
   > P
     V1   V2   V3   V4 
   0.84 0.36 0.36 0.64 

Notice that the SNP locus :math:`k` is most specific to conditions V2 and V3
(0.36), and this is because:

- gene 2 has the lowest specificity percentile (0.2) in condition V2
- gene 4 has the lowest specificity percentile (0.2) in condition V3

