import logging
import os

import pytest

from cellmlmanip import parser

logging.getLogger().setLevel(logging.INFO)


class TestHodgkin:
    @pytest.fixture(scope="class")
    def model(self):
        hodgkin_cellml = os.path.join(
            os.path.dirname(__file__),
            "cellml_files",
            "hodgkin_huxley_squid_axon_model_1952_modified.cellml"
        )
        p = parser.Parser(hodgkin_cellml)
        model = p.parse()
        return model

    def test_counts(self, model):
        # https://models.cellml.org/exposure/5d116522c3b43ccaeb87a1ed10139016/hodgkin_huxley_1952_variant01.cellml/cellml_math
        assert len(model.components) == 8
        eq_count = sum([len(c.equations) for c in model.components.values()])
        assert eq_count == 17

    def test_setup_connections(self, model):
        model.make_connections()
        target = model.components['sodium_channel'].variables['h']
        source = model.components['sodium_channel_h_gate'].variables['h']
        assert target['assignment'] == source['sympy.Dummy']

    def test_add_units_to_equations(self, model):
        for c in model.components.values():
            c.add_units_to_equations()
        equation = model.components['sodium_channel'].equations[0]
        lhs_units = model.units.summarise_units(equation.lhs)
        assert lhs_units == model.units.store['millivolt']

    def test_check_left_right_units(self, model):
        for c in model.components.values():
            for e in c.equations:
                model.check_left_right_units_equal(e)

    def test_get_equations(self, model):
        graph = model.get_equation_graph()

        # use `dot -Tpng path.dot -o path.png`
        # nx.nx_agraph.write_dot(graph,
        #                        '/Users/tamuri/Desktop/path.dot')

        assert len(graph.nodes) == 31

        import networkx as nx
        out = nx.topological_sort(graph)
        for node in out:
            print('%r: %r' % (node, graph.nodes[node]['equation']))
