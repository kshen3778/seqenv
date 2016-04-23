# Built-in modules #

# Internal modules #
from seqenv import module_dir
from seqenv.common.cache import property_cached

# Third party modules #
import networkx, pygraphviz

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
        nodes.update(test_envos)
        nodes = list(nodes)
        # Return #
        return self.networkx.subgraph(nodes)

    def add_weights(self, g, weights=None):
        """Input a networkx DiGraph object.
        Outputs a pygraphviz AGraph object."""
        return networkx.nx_agraph.to_agraph(g)

    def add_style(self, g):
        """Input a pygraphviz AGraph object.
        Outputs a pygraphviz AGraph object."""
        for node in g.nodes():
            node.attr['label'] = node.attr['name']
            node.attr['shape'] = 'box'
            node.attr['style'] = 'rounded'
        for edge in g.edges():
            if edge.attr['label'] == 'located_in': edge.attr['color'] = 'turquoise4'
            edge.attr['label'] = ''
        return g

    def add_legend(self, g):
        """Input a pygraphviz AGraph object.
        Outputs a pygraphviz AGraph object."""
        legend_txt = """
        subgraph legend {
        label = "Legend";
        key [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">
          <tr><td align="right" port="i1">Is a</td></tr>
          <tr><td align="right" port="i2">Part of</td></tr>
          <tr><td align="right" port="i3">Located in</td></tr>
          </table>>];
        key2 [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">
          <tr><td port="i1">&nbsp;</td></tr>
          <tr><td port="i2">&nbsp;</td></tr>
          <tr><td port="i3">&nbsp;</td></tr>
          </table>>];
        key:i1:e -> key2:i1:w [color=red];
        key:i2:e -> key2:i2:w [color=blue];
        key:i3:e -> key2:i3:w [color=turquoise4];
        }"""
        orig_txt = g.to_string().rstrip()[:-1]
        new_text = orig_txt[:-1] + '\n' + legend_txt + '\n' + '}'
        print new_text
        return pygraphviz.AGraph(new_text)

    def write_to_dot(self, g, path):
        """Input a pygraphviz AGraph object."""
        with open(path, 'w') as handle:
            handle.write(g.to_string())

    def draw_to_pdf(self, g, path):
        """Input a pygraphviz AGraph object."""
        g.draw(path, format='pdf', prog='dot')

    # --------------------------- In this section --------------------------- #
    # print_test
    # draw_with_networkx

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
