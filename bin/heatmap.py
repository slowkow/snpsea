#!/usr/bin/env python
"""
:File:     heatmap.py
:Author:   Kamil Slowikowski
:Updated:  July 25, 2013

Create a heatmap of samples and SNPs from the snp_pvalues.txt file produced by
snpsea to reveal the SNPs and genes that contribute to statistically
significant sample specificity scores. Each of the heatmap's cells are
annotated with the given SNP's most specifically expressed gene's name (NCBI
symbol).

Usage:
    heatmap.py --out DIR [--alpha NUM --title STR]

Options:
    -h, --help      Show this message and exit.
    --title STR     Title for the figure.
    --out DIR       Directory with output files created by snpsea.
    --alpha NUM     Significance threshold [default: 0.05].
"""


from docopt import docopt
import matplotlib as mp
import numpy as np
import pandas as pd
import pylab as pl
import textwrap as tw
import re


def main():
    args = docopt(__doc__)
    import os

    out = lambda x: os.path.join(args['--out'], x)
    
    matrix = find_gene_matrix(out('args.txt'))

    heatmap(out('pvalues.txt'),
            matrix,
            out('snp_pvalues.txt'),
            out('snp_pvalues_heatmap.pdf'),
            title=args['--title'] or '',
            alpha=float(args['--alpha']))


def heatmap(f_sample_pvalues, f_matrix, f_pvalues,
            f_heatmap, title='', alpha=0.05):
    sample_pvalues = pd.read_table(f_sample_pvalues)
    # Bonferroni corrected p-values.
    bonf = 1 - (1 - sample_pvalues.pvalue) ** sample_pvalues.shape[0]
    #sample_pvalues.pvalue = bonf
    # Select just the significant samples.
    sample_pvalues = sample_pvalues[bonf < alpha]

    # No sample is significant.
    if sample_pvalues.shape[0] == 0:
        return

    def make_label((sample, pvalue, nulls, reps)):
        if nulls == 0:
            return '{}\np < {}'.format(sample, 1.0 / reps)
        else:
            return '{}\np = {:0.2e}'.format(sample, pvalue)


    samples = sample_pvalues['name'].tolist()[::-1]
    labels = sample_pvalues.apply(make_label, axis=1).tolist()[::-1]

    # Detect compression.
    comp = 'gzip' if f_matrix.endswith('.gz') else None

    # Convert EntrezGene ID numbers to HGNC symbols.
    matrix = pd.read_table(f_matrix, 
                           compression=comp,
                           skiprows=2,
                           index_col=0)
    genelist = list(matrix['Description'].iteritems())
    genedict = {}
    for entrezid, symbol in genelist:
        genedict[entrezid] = symbol

    # Read snp pvalues for each sample and the representative gene symbols.
    df      = pd.read_table(f_pvalues)
    df['gene'] = df['gene'].map(lambda x: genedict.get(x, x))
    # Get a matrix of conditions and SNP loci with gene symbols.
    symbols = df.pivot(index='marker', columns='column', values='gene')
    # And a matrix of conditions and SNP loci with p-values.
    pvalues = df.pivot(index='marker', columns='column', values='pvalue')

    # Select just the significant samples.
    symbols = symbols[samples]
    pvalues = pvalues[samples]

    # Select just SNPs that have at least one significant gene.
    idx = pvalues.apply(lambda row: any(row < alpha), axis=1)
    # No SNP is significant.
    if not any(idx):
        return
    symbols = symbols[idx]
    pvalues = pvalues[idx]
   
    # Sort the SNPs by their max contribution to the significant samples.
    contribs = pvalues.apply(lambda r: np.log(r.replace(0, 1)).min(), axis=1)
    idx = np.argsort(contribs)
    symbols = symbols.ix[idx]
    pvalues = pvalues.ix[idx]

    cell_size = 1.8
    n_snps = pvalues.shape[0]
    n_samples = len(samples)

    fig = pl.figure(figsize=(n_snps * cell_size, 0.6 * n_samples * cell_size))

    # Transform the pvalues for plotting.
    log10_pvalues = -np.log10(pvalues.T)

    minp = np.floor(np.min(log10_pvalues.values))
    maxp = np.ceil(np.max(log10_pvalues.values))
    ticks = np.linspace(minp, maxp, maxp - minp + 1)

    # Plot heatmap.
    ax = fig.add_axes([0.3, 0.1, 0.6, 0.6])
    cmap = pl.cm.Oranges
    norm = mp.colors.BoundaryNorm(ticks, cmap.N)
    im = ax.matshow(log10_pvalues, aspect='auto', origin='lower',
                    cmap=cmap, norm=norm)

    # Ticks.
    ax.tick_params(length=0)
    ax.xaxis.tick_bottom()
    ax.set_xticks(np.arange(n_snps) - 0.45)
    ax.set_yticks(np.arange(n_samples))

    # Labels. If the locus is represented by joined SNPs, print line by line.
    xlabels = [x.replace(',', '\n') for x in list(pvalues.index)]
    ax.set_xticklabels(xlabels, rotation=0, ha='left')
    ax.set_yticklabels(labels)
    ax.set_title(title)

    # Plot the NCBI gene names.
    xy = [(x, y) for x in range(n_snps) for y in range(n_samples)]
    for x, y in xy:
        ax.text(x, y, tw.fill(str(symbols.ix[x, y]), 8),
                bbox=dict(color='white'),
                ha='center', va='center')

    # Plot a colorbar.
    axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.6])

    colorbar = pl.colorbar(im, cax=axcolor, ticks=ticks)

    colorbar.set_label('-log10(p-value)')

    # Save.
    fig.savefig(f_heatmap, bbox_inches='tight')


def find_gene_matrix(filename):
    with open(filename) as f:
        for line in f:
            m = re.search('--gene-matrix\s([^\n]+)', line)
            if m:
                return m.groups()[0].rstrip()
    return None


if __name__ == '__main__':
    main()
