# coding: utf-8

# this file is for those of you with a little more coding experience. You need to calculate
# go-enrichment scores using some gene expression values.

# Feel free to add even more features if you like--the backend code is simple to figure
# out. But don't forget to help your teammates, and to figure out your perturbation!

from collections import Counter, defaultdict

import numpy as np
import scipy.stats as st


import matplotlib.pyplot as plt

# inputs:
# - a list of [(gene name, gene value) ... ]
# - a dictionary of {GOID: [list of genes associated with this term]}
# - an int, N, to use as a parameter
#
# outputs:
# - a list of GOIDs, with associated enrichment scores,
#   testing for positive enrichment (sorted from most significant to least)
# - another list of GOIDs and enrichment scores,
#   testing for negative enrichment (again sorted)
#
def calculate_enrichment(gene_data, go_to_genes, n=100):
    gene_list = sorted(gene_data, key=lambda g: g[1], reverse=True)

    gene_set1 = set(g for g,v in gene_list[:n])
    gene_set2 = set(g for g,v in gene_list[-n:])

    e_scores = dict()
    ne_scores = dict()

    for go in go_to_genes:
        goset = set(go_to_genes[go])
        e_scores[go] = -(st.hypergeom.logsf(len(goset & gene_set1) - 1,
                                            len(gene_list),
                                            len(goset),
                                            len(gene_set1)) / np.log(10))
        ne_scores[go] = -(st.hypergeom.logsf(len(goset & gene_set2) - 1,
                                             len(gene_list),
                                             len(goset),
                                             len(gene_set2)) / np.log(10))

    e_scores = sorted(((go,e_scores[go]) for go in e_scores),
                      key=lambda (go,s): s, reverse=True)
    ne_scores = sorted(((go,ne_scores[go]) for go in ne_scores),
                       key=lambda (go,s): s, reverse=True)

    return e_scores,ne_scores


# You can make your website fancier by creating a figure for each experiment
# (it's up to your what this displays). Install the mpld3 module with pip:
#    pip install mpld3
# and read their website to get started: http://mpld3.github.io/index.html

# uncomment this line to import the module after you've installed it
import mpld3

# input:
# - a list of [(gene name, gene value) ... ]
#
# output:
# - a dictionary, created by the mpld3 module from a figure you've made (see below)
#
def plot_experiment(gene_data):
    # Sort genes by fold expression change
    gene_names,fold_changes = zip(*sorted(gene_data, key=lambda e: e[1], reverse=True))

    # Colors
    almost_black = '#262626'
    ucsf_teal = '#81B0B5'

    # Generate plot in pyplot
    fig, ax = plt.subplots(1)
    scatter = ax.scatter(range(1, len(fold_changes) + 1),
                         fold_changes, label="Datapoint",
                         alpha=0.5, facecolor=ucsf_teal, linewidth=0)

    ax.set_title("Gene Expression Plot", fontweight="bold")
    ax.set_ylabel("Fold Change")
    ax.set_xlabel("Gene Rank")

    ax.set_xlim((-100, len(gene_names) + 100))

    # Pass data to mpld3
    tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=gene_names)
    mpld3.plugins.connect(fig, tooltip)
    mpld3_dict = mpld3.fig_to_dict(fig)

    # Return dict representation of figure
    # This will be converted to JSON in the browser
    return mpld3_dict

