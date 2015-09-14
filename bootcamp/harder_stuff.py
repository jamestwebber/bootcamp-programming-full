# coding: utf-8

# this file is for those of you with a little more coding experience. You need to calculate
# go-enrichment scores using some gene expression values.

# Feel free to add even more features if you like--the backend code is simple to figure
# out. But don't forget to help your teammates, and to figure out your perturbation!

from collections import Counter, defaultdict

import numpy as np
import scipy.stats as st

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

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


# Default distance metrics and clustering methods
DIST_METRICS = set(['correlation', 'cosine', 'euclidean'])
DEF_DIST_METRIC = 'euclidean'
CLUSTER_METHODS = set(['single', 'complete', 'average', 'ward'])
DEF_CLUSTER_METHOD = 'ward'

# dicts to cache gene and experiment clusters using different distance metrics
exp_dist_matrices = {}
gene_dist_matrices = {}
exp_clusters_dict = {}
gene_clusters_dict = {}


# input:
# -  a list or dictionary that maps from the id of an experiment (an int: 0, 1, ..)
#    to a list of (systematic name, fold-change value) tuples
# outputs:
# - NxM numpy array with N experiments and M genes
# - list of experiment numbers
# - list of gene names
def make_exp_array(exp_data):
    values = []
    exp_labels = []
    gene_labels = []

    for exp_num in range(len(exp_data)):
        _labels, _values = zip(*exp_data[exp_num])
        values.append(_values)
        exp_labels.append(exp_num)
    gene_labels = _labels

    return np.asarray(values, dtype=np.float), exp_labels, gene_labels


# input:
# -  an NxM data array
# outputs:
# - square distance matrix along axis 0
def calculate_distance_matrix(data_array, metric=DEF_DIST_METRIC):
    return squareform(pdist(data_array, metric))


# input:
# -  a list or dictionary that maps from the id of an experiment (an int: 0, 1, ..)
#    to a list of (systematic name, fold-change value) tuples
# output:
# - a dictionary, created by the mpld3 module from a figure you've made (see below)
def plot_experiment_clusters(exp_data, metric=DEF_DIST_METRIC,
                             method=DEF_CLUSTER_METHOD):
    global exp_dist_matrices
    global gene_dist_matrices
    global exp_clusters_dict
    global gene_clusters_dict
    if metric not in DIST_METRICS:
        metric = DEF_DIST_METRIC
    if method not in CLUSTER_METHODS:
        method = DEF_CLUSTER_METHOD

    print("Making data array")
    exp_data_array, exp_labels, gene_labels = make_exp_array(exp_data)
    if metric not in exp_dist_matrices:
        print("Making exp dist matrix")
        exp_dist_matrices[metric] = calculate_distance_matrix(
            exp_data_array)
    if (metric, method) not in exp_clusters_dict:
        print("Making exp clusters")
        exp_clusters_dict[(metric, method)] = linkage(
            exp_dist_matrices[metric], method=method)
    if metric not in gene_dist_matrices:
        print("Making gene dist matrix")
        gene_dist_matrices[metric] = calculate_distance_matrix(
            exp_data_array.T)
    if (metric, method) not in gene_clusters_dict:
        print("Making gene clusters")
        gene_clusters_dict[(metric, method)] = linkage(
            gene_dist_matrices[metric], method=method)
    print("Plotting clusters")

    # Compute and plot first dendrogram.
    fig = plt.figure(figsize=(8, 8))

    # Compute and plot experiments dendrogram.
    ax1 = fig.add_axes([0.3, .70, 0.6, 0.30])
    exp_dendro = dendrogram(exp_clusters_dict[(metric, method)],
                            labels=exp_labels)
    idx1 = exp_dendro['leaves']
    ax1.set_xticks([])
    ax1.set_yticks([])

    # Compute and plot genes dendrogram.
    ax2 = fig.add_axes([0, 0.09, 0.30, 0.61])
    gene_dendro = dendrogram(gene_clusters_dict[(metric, method)],
                             labels=gene_labels, orientation='right')
    idx2 = gene_dendro['leaves']
    ax2.set_xticks([])
    ax2.set_yticks([])

    # Plot heatmap of experimental results
    axmatrix = fig.add_axes([0.30, 0.09, 0.6, 0.61])
    sorted_data_matrix = np.asarray([[exp_data_array[x][y] for y in idx2]
                                     for x in idx1],
                                    dtype=np.float).T

    im = axmatrix.matshow(sorted_data_matrix, aspect='auto', origin='lower',
                          cmap="PuBuGn")
    sorted_exp_labels = [exp_labels[x] for x in idx1]
    axmatrix.set_xticklabels(sorted_exp_labels)
    axmatrix.tick_params(axis="both", which="both", bottom="off", top="off",
                         labelbottom="on", left="off", right="off",
                         labelleft="off")
    axmatrix.set_yticks([])
    axmatrix
    axcolor = fig.add_axes([0.91, 0.09, 0.02, 0.61])
    fig.colorbar(im, cax=axcolor)

    mpld3_dict = mpld3.fig_to_dict(fig)

    return mpld3_dict
