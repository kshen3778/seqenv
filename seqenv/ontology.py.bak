# Built-in modules #

# Internal modules #
from seqenv import module_dir
from seqenv.common.cache import property_cached

# Third party modules #
import sh, networkx
import matplotlib.colors

# A list of envos to help test this module #
test_envos = [
    "ENVO:00000033",
    "ENVO:00000043",
    "ENVO:00000067",
    "ENVO:00000143",
    "ENVO:00000210",
    "ENVO:00000215",
    "ENVO:00000475",
]

################################################################################
class Ontology(object):
    """A object that gives you access to the graph (network with nodes and edges)
    of the ENVO ontology from the OBO file's path.

    Other libraries not used here that could be added:
        * graphviz:   http://graphviz.readthedocs.org/en/latest/api.html#digraph
        * pydot:      https://github.com/erocarrera/pydot
    """

    def __init__(self, path=None):
        """Give the path to the OBO file"""
        if path is None: path = module_dir + 'data_envo/envo.obo'
        self.path = path

    # --------------------------- In this section --------------------------- #
    # orange_obo
    # goatools
    # orange_network
    # pygraphviz
    # networkx

    @property_cached
    def orange_obo(self):
        """The ontology loaded by the `orange` library.
        * http://orange.biolab.si
        * http://orange-bioinformatics.readthedocs.org/en/latest/
        * https://github.com/biolab/orange-bio
        * https://bitbucket.org/biolab/orange-bioinformatics
        To install: $ pip install Orange-Bioinformatics
        """
        from orangecontrib.bio.ontology import OBOOntology
        return OBOOntology(self.path)

    @property_cached
    def goatools(self):
        """The network loaded into goatools' format.
        * https://github.com/tanghaibao/goatools
        To install: $ pip install goatools
        """
        from goatools import obo_parser
        return obo_parser.GODag(self.path)

    @property_cached
    def orange_network(self):
        """The network converted to `orange network` format.
        Doesn't seem to work until they update PyPI.
        * https://bitbucket.org/biolab/orange-network/
        * http://orange-network.readthedocs.org/en/latest/
        To install: $ pip install orange-network
        """
        return self.orange_obo.to_network()

    @property_cached
    def pygraphviz(self):
        """The network converted to `pygraphviz` format.
        * http://pygraphviz.github.io/documentation/pygraphviz-1.3rc1/
        To install: $ pip install pygraphviz
        """
        g = self.orange_obo.to_graphviz()
        assert g.is_directed()
        assert g.is_strict()
        return g

    @property_cached
    def networkx(self):
        """The network converted to `networkx` format.
        Seems like it looses directionality.
        * https://networkx.readthedocs.org/en/stable/
        To install: $ pip install networkx
        """
        g = self.orange_obo.to_networkx()
        assert networkx.is_directed_acyclic_graph(g)
        return g

    # --------------------------- In this section --------------------------- #
    # test
    # get_subgraph
    # add_weights
    # draw_to_pdf
    # write_to_dot

    def get_subgraph(self, envos=None):
        """Given a list of ENVO terms, get the subgraph that contains them all
        and all their ancestors, up to the root.
        Outputs a networkx DiGraph object."""
        # Testing mode #
        if envos is None: envos = test_envos
        # All nodes #
        nodes = set(n for e in envos for n in networkx.descendants(self.networkx, e))
        nodes.update(envos)
        nodes = list(nodes)
        # Return #
        return self.networkx.subgraph(nodes)

    def add_weights(self, g, weights=None):
        """Input a networkx DiGraph object.
        Outputs a pygraphviz AGraph object."""
        g = networkx.nx_agraph.to_agraph(g)
        if weights is None: return g
        for envo in weights:
            node   = g.get_node(envo)
            weight = weights[envo]
            color  = matplotlib.colors.rgb2hex((1.0, 1.0 - weight, 0.0))
            node.attr['fillcolor'] = color
        return g

    def add_style(self, g):
        """Input a pygraphviz AGraph object.
        Outputs a pygraphviz AGraph object."""
        for node in g.nodes():
            text = node.attr['name']
            node.attr['label'] = text.replace(' ','\\n')
            node.attr['name']  = ''
            node.attr['shape'] = 'Mrecord'
            node.attr['style'] = 'filled'
            # To add the envo id to each node, uncomment:
            #envo = node.attr['label']
            #node.attr['label'] = "{<f0> %s|<f1> %s}" % (envo, text)
        for edge in g.edges():
            if edge.attr['label'] == 'located_in': edge.attr['color'] = 'turquoise4'
            edge.attr['label'] = ''
        return g

    def write_to_dot(self, g, path):
        """Input a pygraphviz AGraph object."""
        with open(path, 'w') as handle: handle.write(g.to_string())

    def add_legend(self, path):
        """Input the path to a dot file."""
        legend_txt = """
        digraph {
          rankdir=LR
          node [shape=plaintext,fontname="helvetica"]
          subgraph cluster_01 {
          label = "NB: darker nodes weigh more";
            key [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">
              <tr><td align="right" port="i1">Is</td></tr>
              <tr><td align="right" port="i2">Part</td></tr>
              <tr><td align="right" port="i3">Located</td></tr>
              </table>>];
            key2 [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">
              <tr><td port="i1">a</td></tr>
              <tr><td port="i2">of</td></tr>
              <tr><td port="i3">in</td></tr>
              </table>>];
            key:i1:e -> key2:i1:w [color=red];
            key:i2:e -> key2:i2:w [color=blue];
            key:i3:e -> key2:i3:w [color=turquoise4];
        }"""
        orig_txt = [line.rstrip('\n') for line in open(path, 'r') if line]
        new_text = [line.lstrip() for line in legend_txt.split('\n') if line]
        new_text = '\n'.join(new_text + orig_txt[2:])
        with open(path, 'w') as handle: handle.write(new_text)

    def draw_to_pdf(self, in_path, out_path):
        """Input a path to a dot file."""
        sh.dot(in_path, '-Tpdf', '-o', out_path)

    # --------------------------- In this section --------------------------- #
    # descends

    def descends(self, e, root):
        """Does the envo term `e` descend from the node `root`?
        Returns True or False."""
        # Auto conversion #
        if isinstance(e,    int):    e = "ENVO:%08d" % e
        if isinstance(root, int): root = "ENVO:%08d" % root
        # Return #
        return e in networkx.ancestors(self.networkx, root)

    # --------------------------- In this section --------------------------- #
    # print_test
    # draw_with_networkx
    # draw_with_pygraphviz

    def print_test(self, e=None):
        """Just a method to see a bit how the different libraries work."""
        # Test node #
        if e is None: e = test_envos[0]
        # Goa #
        print "Goa: "
        print self.goatools[e]
        # Pygraphviz #
        print "pygraphviz: "
        print self.pygraphviz[e]
        print self.pygraphviz.successors(e)
        print self.pygraphviz.predecessors(e)
        print self.pygraphviz.get_node(e)
        # Networkx #
        import networkx
        print "networkx: "
        print self.networkx[e]
        print self.networkx.successors(e)
        print self.networkx.predecessors(e)
        print networkx.ancestors(self.networkx, e)   # same as predecessors
        print networkx.descendants(self.networkx, e) # almost as child_to_parents

    def draw_with_networkx(self, g, path):
        """Input a networkx DiGraph object."""
        from matplotlib import pyplot
        networkx.draw(g)
        pyplot.savefig(path)
        pyplot.close()

    def draw_with_pygraphviz(self, g, path):
        """Input a pygraphviz AGraph object."""
        with open(path, 'w') as handle:
            handle.write(g.to_string())

