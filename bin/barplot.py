#!/usr/bin/env python
"""
:File:     barplot.py
:Author:   Kamil Slowikowski
:Updated:  September 21, 2013

Create a bar plot of p-values for each sample with an adjacent heatmap
triangle showing all pairwise sample Spearman correlations of expression
profiles.

Usage:
    barplot.py --out DIR [--top INT --title STRING --width INT --fontsize INT]

Options:
    -h, --help         Show this message and exit.
    --out DIR          Directory with SNPspec output files.
    --top INT          Only plot the top N results [default: 25].
    --title STRING     Title printed at the top of the figure.
    --width INT        Width of the figure [default: 12].
    --fontsize INT     Font size of y-axis labels [default: 12].
"""


from docopt import docopt
import itertools
import matplotlib as mp
import numpy as np
import pandas as pd
import pylab as pl
import scipy.cluster.hierarchy as sch
import os
import re


def main():
    args = docopt(__doc__)
    import os

    out = lambda x: os.path.join(args['--out'], x)
   
    width = int(args['--width'])
    height = width / 2

    matrix = find_gene_matrix(out('args.txt'))

    if not args['--title']:
        args['--title'] = os.path.basename(args['--out'])

    barplot(out('pvalues.txt'),
            matrix,
            out('snp_genes.txt'),
            out('pvalues_barplot.pdf'),
            title=args['--title'],
            figsize=(width, height),
            fontsize=int(args['--fontsize']),
            top=int(args['--top']))


def barplot(f_pvalues, f_matrix, f_genes, f_plot, figsize=(5, 5), fontsize=12,
            cluster_method='all', title=None, alpha=0.05, top=50):
    # Read results from SNPspec.
    pvalues = pd.read_table(f_pvalues, index_col=0)

    top = min(top, pvalues.shape[0])

    # Get rid of the 0 p-values.
    idx = pvalues.pvalue == 0
    pvalues.ix[idx, 'pvalue'] = 1.0 / pvalues.ix[idx, 'nulls_tested']
    pvalues = pvalues.sort('pvalue')[:top]

    # Detect compression.
    comp = 'gzip' if f_matrix.endswith('.gz') else None

    # The matrix must be in GCT format.
    matrix = pd.read_table(f_matrix, 
                               compression=comp,
                               skiprows=2,
                               index_col=0).drop('Description', axis=1)

    # Select the top N columns.
    n_samples = matrix.shape[1]
    matrix = matrix.ix[:, pvalues.index]

    # By default, cluster samples by all of the gene values.
    # Otherwise, cluster by the genes overlapped by the user's SNPs.
    if cluster_method == 'user':
        snp_genelists = pd.read_table(f_genes)
        splitcomma = lambda x: str(x).split(',')
        genelists = snp_genelists['genes'].dropna().apply(splitcomma)
        snp_genes = map(int, set(flatten(genelists)))
        matrix = matrix.ix[snp_genes, :]

    # Create a figure.
    fig = pl.figure(figsize=figsize)

    # Axes for the heatmap triangle.
    ax = fig.add_subplot(121, frame_on=False, aspect=2.0)

    # Get the heatmap triangle's axes and the order of the clustered samples.
    cax, order = heatmap_triangle(matrix, ax)

    # Adjust spacing between the heatmap triangle and the barplot.
    fig.subplots_adjust(wspace=0, hspace=0, left=0, right=0.4)

    # Axes for the barplot.
    ax = fig.add_subplot(122, frame_on=False)

    # Order the p-values by the clustering.
    pvalues = pvalues.ix[order]
    
    # Magic fontsize formula obtained by trial and error.
    #fontsize = 7 * figsize[1] ** (figsize[1] / float(top))

    # Plot the negative log10 transformed p-values.
    pvalues['pvalue'] = -np.log10(pvalues['pvalue'])
    ax = pvalues['pvalue'].plot(ax=ax,
                                kind='barh',
                                title=title,
                                linewidth=0,
                                #fontsize=fontsize,
                                grid=False,
                                color='blue')

    # Vertical grid lines.
    ax.grid(b=True, which='major', axis='both', alpha=0.5)

    # Tick marks for the x-axis.
    xticks = np.arange(0, round(max(pvalues['pvalue'])) + 1)
    ax.set_xticks(xticks)

    # Put the y-axis marks on the right.
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')

    # Adjust tick parameters.
    ax.tick_params(length=0, axis='x')
    ax.tick_params(length=0, axis='y', labelsize=fontsize)

    # Labels.
    ax.set_xlabel('-log10(p-value)')
    ax.set_ylabel('')

    # Add a vertical Bonferroni p-value cutoff line.
    bonf = lambda a, n: 1 - (1 - a) ** (1 / float(n))

    ax.vlines(x=-np.log10(bonf(alpha, n_samples)),
              ymin=min(ax.get_yticks()) - 0.5,
              ymax=max(ax.get_yticks()) + 0.5,
              color='red')

    # Save.
    fig.savefig(f_plot, bbox_inches='tight')


def heatmap_triangle(dataframe, axes):
    """Create a heatmap of the lower triangle of a pairwise correlation
    matrix of all pairs of columns in the given dataframe. The heatmap
    triangle is rotated 45 degrees clockwise and drawn on the given axes.

    :param dataframe:   pandas.DataFrame
    :param axes:        matplotlib.axes.Axes
    """
    N = dataframe.shape[1]
    D = dataframe.corr(method='pearson')

    # UPGMA clustering, but other methods are also available.
    Z = sch.linkage(D, method='average')
    R = sch.dendrogram(Z, no_plot=True)
    sample_order = R['leaves']
    D = D.ix[sample_order, sample_order]

    # Get the lower triangle of the matrix. 
    C = np.tril(D)#, -1)
    # Mask the upper triangle and the diagonal.
    C = np.ma.masked_array(C, C == 0)
    for i in range(N):
        C[i, i] = 0

    # Transformation matrix for rotating the heatmap.
    t = np.array([[0.5, 1], [0.5, -1]])

    B = itertools.product(range(N, -1, -1), range(0, N + 1, 1))
    A = np.array([(i[1], i[0]) for i in B])
    A = np.dot(A, t)

    # -1.0 correlation is blue, 0.0 is white, 1.0 is red.
    cmap = pl.cm.RdBu_r
    norm = mp.colors.BoundaryNorm(np.linspace(-1, 1, 14), cmap.N)

    # This MUST be before the call to pl.pcolormesh() to align properly.
    axes.set_xticks([])
    axes.set_yticks([])

    # Plot the correlation heatmap triangle.
    X = A[:, 1].reshape(N + 1, N + 1)
    Y = A[:, 0].reshape(N + 1, N + 1)
    caxes = pl.pcolormesh(X, Y, np.flipud(C), axes=axes, cmap=cmap, norm=norm)

    # Remove the ticks and reset the x limit.
    axes.set_xlim(right=0)

    # Add a colorbar below the heatmap triangle.
    pl.colorbar(caxes, ax=axes, orientation='horizontal', shrink=0.7,
                fraction=0.05, pad=-0.05, ticks=np.linspace(-1, 1, 5))

    return caxes, D.index


def flatten(lists):
    """Flatten a list of lists into a single list.

    :param lists:  a list

    >>> flatten([1, 2, [3, 4], 5, [6, [7, [8, 9]]]])
    [1, 2, 3, 4, 5, 6, 7, 8, 9]
    """
    return list(itertools.chain.from_iterable(lists))


def find_gene_matrix(filename):
    with open(filename) as f:
        for line in f:
            m = re.search('--gene-matrix\s([^\n]+)', line)
            if m:
                return m.groups()[0].rstrip()
    return None


if __name__ == '__main__':
    main()
