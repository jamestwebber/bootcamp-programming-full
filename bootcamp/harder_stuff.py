# coding: utf-8

# this file is for those of you with a little more coding experience. You need to calculate
# go-enrichment scores using some gene expression values.

# Feel free to add even more features if you like--the backend code is simple to figure
# out. But don't forget to help your teammates, and to figure out your perturbation!

import easier_stuff as es

from collections import Counter, defaultdict

import numpy as np
import scipy.stats as st

import scipy.cluster.hierarchy as sch

import matplotlib.pyplot as plt

import matplotlib.figure
import matplotlib.colors
from matplotlib.gridspec import GridSpec

import easier_stuff as es

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


# inputs:
# - the index for a single experiment (an int: 0, 1, ..)
# - a list or dictionary that maps from id to experiment data--this will be
#   the output from easier_stuff.experiment()
# - an int, N, to use as a parameter
#
# outputs:
# - a list of the N most similar experiments (determined however you like)
#
def similar_experiments(exp, exp_data, n=10):
    # this assumes data for every experiment is in the same order...
    comp_e = lambda i: 1.0 - st.pearsonr([v for g,v in exp_data[exp]],
                                         [v for g,v in exp_data[i]])[0]

    most_similar = sorted([i for i in range(len(exp_data)) if i != exp],
                          key=comp_e)[:n]

    return most_similar


def go2gene_network(goid, n, graph=None):
    if graph is None:
        graph = {'nodes': [{'id': goid, 'node_type': es.go_info(goid)[1], 'name': es.go_info(goid)[0]}],
                 'links': [],
                 'nodeset': {goid}}

    node_i = len(graph['nodes']) - 1

    if n <= 0:
        return graph

    gene_list = es.go_to_gene(goid)

    if len(gene_list) > 100:
        graph['links'].append({'source': node_i, 'target': len(graph['nodes'])})
        graph['nodes'].append({'id': '', 'node_type': 'gene',
                               'name': '({} genes)'.format(len(gene_list))})
        return graph

    for gene in gene_list:
        if gene in graph['nodeset']:
            continue
        graph['links'].append({'source': node_i, 'target': len(graph['nodes'])})
        graph['nodes'].append({'id': gene, 'node_type': 'gene', 'name': es.gene_name(gene)})
        graph['nodeset'].add(gene)
        gene2go_network(gene, n - 1, graph)

    return graph


def gene2go_network(gene, n, graph=None):
    if graph is None:
        graph = {'nodes': [{'id': gene, 'node_type': 'gene', 'name': es.gene_name(gene)}],
                 'links': [],
                 'nodeset': {gene}}

    node_i = len(graph['nodes']) - 1

    if n <= 0:
        return graph

    for goid in es.gene_to_go(gene):
        if goid in graph['nodeset']:
            continue
        graph['nodes'].append({'id': goid, 'node_type': es.go_info(goid)[1], 'name': es.go_info(goid)[0]})
        graph['links'].append({'source': node_i, 'target': len(graph['nodes']) - 1})
        graph['nodeset'].add(goid)
        go2gene_network(goid, n - 1, graph)

    return graph


# inputs:
# - a goid or gene systematic name
# - a flag that is True when the previous input is a goid
# - an int n that describes how many steps in the network
#
# note: your browser will get cranky if there are too many nodes
#
# outputs:
# - a dictionary containing two keys: 'nodes' and 'links'
#   'nodes' : a list of node dictionaries, of the form:
#             {'id': goid or gene systematic name,
#               'name': go or gene name
#               'node_type': go aspect or 'gene'}
#   'links' : a list of edge dictionaries, of the form:
#             {'source': index of the source node,
#              'target': index of the target node}
#   where the indexes point to nodes in the 'nodes' list
#
#  e.g. go_network('GO:0006383') returns
#  { 'nodes': [{'id': u'GO:0006383',
#               'name': 'transcription from RNA polymerase III promoter',
#               'node_type': 'P'},
#              {'id': 'YAL001C', 'name': 'TFC3', 'node_type': 'gene'},
#              {'id': 'YBR123C', 'name': 'TFC1', 'node_type': 'gene'},
#              ...
#             ]
#    'links': [{'source': 0, 'target': 1},
#              {'source': 0, 'target': 2},
#              {'source': 0, 'target': 3}
#              ...
#             ]
def go_network(goid_or_gene, is_goid=True, n=2):
    if is_goid:
        gnet = go2gene_network(goid_or_gene, n)
    else:
        gnet = gene2go_network(goid_or_gene, n)

    return {'nodes': gnet['nodes'], 'links': gnet['links']}


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

import os
import pickle

if os.path.exists('clustering.pickle'):
    with open('clustering.pickle') as f:
        _exp_matrix,Y1,Y2,idx1,idx2,_gene_labels = pickle.load(f)
    print "loaded pickle"
else:

    _gene_labels = [g for g,v in es.experiment()[0]]
    _exp_matrix = np.vstack([[v for g,v in es.experiment()[i]]
                             for i in sorted(es.experiment())])
    print "clustering experiments"
    Y1 = sch.linkage(_exp_matrix, method='ward', metric='euclidean')
    print "clustering genes"
    Y2 = sch.linkage(_exp_matrix.T, method='ward', metric='euclidean')
    print "done"
    idx1, idx2 = sch.leaves_list(Y1), sch.leaves_list(Y2)
    _gene_labels = [_gene_labels[i] for i in idx2]
    _exp_matrix = _exp_matrix[np.ix_(idx1, idx2)]

    with open('clustering.pickle', 'w') as OUT:
        pickle.dump((_exp_matrix, Y1, Y2, idx1, idx2, _gene_labels), OUT)

# input:
# -  a list or dictionary that maps from the id of an experiment (an int: 0, 1, ..)
#    to a list of (systematic name, fold-change value) tuples
# output:
# - a dictionary, created by the mpld3 module from a figure you've made (see below)
def plot_experiment_clusters():
    fig = plt.figure(figsize=(8.,8.))

    gs = GridSpec(2, 3, left=0.05, bottom=0.05, right=0.95, top=0.95,
                  wspace=0.05, hspace=0.05, height_ratios=(1,2), width_ratios=(4,8,1))

    ax01 = fig.add_subplot(gs[1])
    ax01.grid('off')
    ax10 = fig.add_subplot(gs[3])
    ax10.grid('off')

    sch.dendrogram(Y2, ax=ax01)
    ax01.set_axis_off()
    ax01.tick_params(bottom="off", labelbottom='off',
                     top="off", labeltop='off',
                     left="off", labelleft='off',
                     right="off", labelright='off')

    sch.dendrogram(Y1, orientation='right', ax=ax10)
    ax10.set_axis_off()
    ax10.tick_params(bottom="off", labelbottom='off',
                     top="off", labeltop='off',
                     left="off", labelleft='off',
                     right="off", labelright='off')

    ax11 = fig.add_subplot(gs[4])

    vmin,vmax = np.percentile(_exp_matrix, (1.0, 99.0))
    norm = matplotlib.colors.Normalize(vmin, vmax)

    im = ax11.matshow(_exp_matrix, aspect='auto', origin='lower',
                      norm=norm, cmap='seismic')

    ax11.tick_params(bottom="off", labelbottom='off',
                     top="off", labeltop='off',
                     left="off", labelleft='off',
                     right="off", labelright='off')

    ax12 = fig.add_subplot(gs[5])
    fig.colorbar(im, cax=ax12)
    ax12.tick_params(bottom='off')
    ax12.set_yticklabels([vmin, vmax])

    mpld3_dict = mpld3.fig_to_dict(fig)

    return mpld3_dict
