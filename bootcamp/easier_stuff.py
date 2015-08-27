# coding: utf-8

# This file is for those of you who learned Python over the summer (you did that, right?).
# In this file, I've put all of the nitty-gritty details of what makes this website work.

# Except it doesn't work, because you need to write all the functions!

# Some of these functions will just make the website easier to use. Some of them are
# important for the enrichment and clustering tasks that your teammates are working on.

# If you need any help, ask your team or a TA.


# (don't delete this but don't worry about it either)
import os # a built-in module, for dealing with filenames
from . import app # this is part of the website guts



# These are all the files you have to work with. Go open them in a text  editor so you can
# get a feel for what they look like, because you need to parse each one to turn on a
# piece of the website.

# A list of yeast genes, with standard names and short descriptions.
GENE_INFO = os.path.join(app.root_path, 'data', 'gene_info.txt')

# A file that maps from GOID to name, aspect (process/function/component), etc
GO_INFO = os.path.join(app.root_path, 'data', 'go_info.txt')

# A two-column file that maps GOID to yeast genes
GO_MEMBERSHIP = os.path.join(app.root_path, 'data', 'go_membership.txt')

# A many-columned file that contains experimental data (yeast microarrays). Each column
# (after the first) is a different experiment, and each row is a gene. The values are log2
# ratios versus untreated control.
EXPERIMENT_FILE = os.path.join(app.root_path, 'data', 'experiment_data.txt')

import csv
from collections import defaultdict

with open(EXPERIMENT_FILE) as f:
    rdr = csv.reader(f, delimiter='\t')
    h = rdr.next()[1:]
    rows = list(rdr)

    _exp_data = dict()
    _exp_data_by_gene = defaultdict(list)

    for i,exp_name in enumerate(h):
        _exp_data[i] = [(row[0],float(row[i+1])) for row in rows]

        for g,v in _exp_data[i]:
            _exp_data_by_gene[g].append(v)

with open(GENE_INFO) as f:
    rdr = csv.DictReader(f, delimiter='\t')

    _gene_names = dict()
    _gene_info = dict()

    for row in rdr:
        _gene_names[row['systematic_name']] = row['gene_name']
        _gene_info[row['systematic_name']] = row['gene_description']


with open(GO_INFO) as f:
    rdr = csv.DictReader(f, delimiter='\t')

    _go_info = dict()
    _go_aspect = defaultdict(set)

    for row in rdr:
        _go_info[row['goid']] = (row['go_term'], row['go_aspect'], row['go_definition'])
        _go_aspect[row['go_aspect']].add(row['goid'])

    _go_aspect = {g:sorted(_go_aspect[g]) for g in _go_aspect}


with open(GO_MEMBERSHIP) as f:
    rdr = csv.DictReader(f, delimiter='\t')

    _gene_to_go = defaultdict(set)
    _go_to_gene = defaultdict(set)

    for row in rdr:
        _gene_to_go[row['systematic_name']].add(row['goid'])
        _go_to_gene[row['goid']].add(row['systematic_name'])

    _gene_to_go = {g:sorted(_gene_to_go[g]) for g in _gene_to_go}
    _go_to_gene = {g:sorted(_go_to_gene[g]) for g in _go_to_gene}


# return a list or dictionary that maps from the id of an experiment (an int: 0, 1, ..)
# to a list of (systematic name, fold-change value) tuples
# e.g. [[('YAL001C', -0.06), ('YAL002W', -0.3), ('YAL003W', -0.07), ... ],
#       [('YAL001C', -0.58), ('YAL002W', 0.23), ('YAL003W', -0.25), ... ],
#        ... ]
def experiment():
    return _exp_data


# map from a gene's systematic name to its standard name
# e.g. gene_name('YGR188C') returns 'BUB1'
def gene_name(gene):
    return _gene_names[gene]


def gene_data(gene):
    return _exp_data_by_gene[gene]


# map from a systematic name to some info about the gene (whatever you want),
# e.g  'YGR188C' -> 'Protein kinase involved in the cell cycle checkpoint into anaphase'
def gene_info(gene):
    return _gene_info[gene]


# map from a systematic name to a list of GOIDs that the gene is associated with
# e.g. 'YGR188C' -> ['GO:0005694', 'GO:0000775', 'GO:0000778', ... ]
def gene_to_go(gene):
    return _gene_to_go[gene]


# map from one of the GO aspects (P, F, and C, for Process, Function, Component),
# to a list of all the GOIDs in that aspect
# e.g. 'C' -> ['GO:0005737', 'GO:0005761', 'GO:0005763', ... ]
def go_aspect(aspect):
    return _go_aspect[aspect]


# map from a GOID (e.g. GO:0005737) to a *tuple* of the term and term definition
# e.g. 'GO:0005737' -> ('cytoplasm', 'All of the contents of a cell... (etc)'
def go_info(goid):
    return _go_info[goid]


# the reverse of the gene_to_go function: map from a GOID
# to a list of genes (systematic names)
# e.g. 'GO:0005737' -> ['YAL001C', 'YAL002W', 'YAL003W', ... ]
def go_to_gene(goid):
    return _go_to_gene[goid]
