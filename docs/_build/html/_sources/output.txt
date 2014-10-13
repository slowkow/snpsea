Output Visualizations
---------------------

View enrichment of tissue-specific gene expression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: figures/Red_blood_cell_count-Harst2012-45_SNPs-GeneAtlas2004-single-pvalues_barplot.png
   :alt: 

A horizontal bar plot of negative log10 p-values for a test of 37
LDL-associated SNPs for enrichment of tissue-specific expression in
profiles of 79 human tissues.

.. code:: bash

    python bin/snpsea-barplot out

View the most specifically expressed gene for each SNP-tissue pair
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: figures/Red_blood_cell_count-Harst2012-45_SNPs-GeneAtlas2004-single-snp_condition_heatmap.png
   :alt: 

A heatmap exposing the contributions of specifically expressed genes
within each SNP linkage interval to the specificity scores of each
tissue.

.. code:: bash

    python bin/snpsea-heatmap out

View the type 1 error rate estimates for each tissue
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: figures/type1error_GeneAtlas2004.png
   :alt: 

A scatter plot of the observed proportion of p-values under various
thresholds after repeating the analysis with 10,000 random SNP sets.

.. code:: bash

    Rscript bin/snpsea-type1error out

