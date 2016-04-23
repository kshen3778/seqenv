"""
# How it worked before #

The scripts from a file that is a transposed csv version of the current output file
'seq_to_concepts.tsv' an example of the transformed is included here 'stc.csv'
The first of the old scripts splits this into seperate files one for each sequence

-> perl ./envo_extract_records.pl -i stc.csv -l sequences_list.txt

* The file sequences_list.txt stores the list of file names
 Then the dot file generating script is run for each sequence file

-> perl ./envo_gen_dot_files.pl -i "$id" -o envo.obo 2>/dev/null

Each dot file as I understand comprises that subset of the obo graph that is linked from an obo root to each of the concepts. With a node whose color reflects the importance of the concept.

# Idea: Try with the parents relationship file.

# Potential client views #
* Grpahviz:
* Gephi:
* Tulip:

# Pip installation #
 pip install graphviz
 pip install pydot
 pip install pygraphviz
 pip install networkx
 pip install Orange-Bioinformatics
 pip install orange-network

# Implementation #
 OK so it seems to me it's like this:
 Given the directed network "n" (the whole onotology)
 Given the three nodes (envo terms) of interest for sequence at hand S which are E1, E2, E3
 Given the four roots which are: R1, R2, R3, R4
 the result is:
 r = n.all_forward_paths(E1, R1) + n.all_forward_paths(E1, R2) + n.all_forward_paths(E1, R3) + n.all_forward_paths(E1, R4) +
 n.all_forward_paths(E2, R1) + n.all_forward_paths(E2, R2) + n.all_forward_paths(E2, R3) + n.all_forward_paths(E2, R4) +
 n.all_forward_paths(E3, R1) + n.all_forward_paths(E3, R2) + n.all_forward_paths(E3, R3) + n.all_forward_paths(E3, R4)
And then r.colorize(E1, E1.weight) etc.

# Graphviz graph:
-> graphviz:   http://graphviz.readthedocs.org/en/latest/api.html#digraph
-> pydot:      https://github.com/erocarrera/pydot
-> pygraphviz: http://pygraphviz.github.io/documentation/pygraphviz-1.3rc1/

* Networkx graph: https://networkx.readthedocs.org/en/stable/
* Orange   graph: http://orange-network.readthedocs.org/en/latest/
->orange-network2 https://bitbucket.org/biolab/orange-network/
->orange-network3 https://bitbucket.org/biolab/orange-network/

* Goa: https://github.com/tanghaibao/goatools
"""

# Modules #
import networkx, os

# Paths #
home = os.environ.get('HOME', '~') + '/'
obo_path = home + 'repos/seqenv/seqenv/data_envo/envo.obo'

# Object #
from seqenv.ontology import Ontology
ontology = Ontology()

# Node test #
from seqenv.ontology import test_envos
nodes = set(n for e in test_envos for n in ontology.networkx.successors(e))
nodes.update(test_envos)
nodes = list(nodes)
subgraph = ontology.networkx.subgraph(nodes)

# Draw test #
ontology.draw_with_networkx(subgraph, home+'Desktop/viz/subgraph.pdf')
converted = networkx.nx_agraph.to_agraph(subgraph)
converted.draw(home+'Desktop/viz/final.pdf', format='pdf', prog='dot')
ontology.write_to_dot(converted, home+'Desktop/viz/final.dot')

# Complete test #
ontology.draw_with_networkx(ontology.networkx, home+'Desktop/viz/networkx_complete.pdf')
ontology.pygraphviz.draw(home+'Desktop/viz/graphviz_complete.pdf', prog='dot')
converted = networkx.nx_agraph.to_agraph(ontology.networkx)
converted.draw(home+'Desktop/viz/converted_complete.pdf', prog='dot')

################################################################################
# With class #
import os
home = os.environ.get('HOME', '~') + '/'
from seqenv.ontology import Ontology
ontology = Ontology()
g = ontology.get_subgraph()
g = ontology.add_weights(g)
g = ontology.add_style(g)
ontology.write_to_dot(g, home+'Desktop/viz/final.dot')
ontology.add_legend(home+'Desktop/viz/final.dot')
ontology.draw_to_pdf(home+'Desktop/viz/final.dot', home+'Desktop/viz/final.pdf')

################################################################################
# Old method #
import pygraphviz
from seqenv import onotology
def build_with_parents(self, envos=None):
    """Given a list of ENVO terms, build a simple network with
    all the keys and all the parents."""
    # Testing mode #
    if envos is None: envos = onotology.test_envos
    # New graph #
    g = pygraphviz.AGraph(directed=True)
    # Main loop #
    for e in envos:
        g.add_node(e)
    # Return #
    return g

