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
        eq_count = len(list(model.equations()))
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
        for e in model.equations():
            model.check_left_right_units_equal(e)

    def test_get_equations(self, model):
        graph = model.get_equation_graph()
        assert len(graph.nodes) == 32

        free_variable = model.find_variable({'type': 'free'})
        assert len(free_variable) == 1
        assert free_variable[0]['cmeta:id'] == 'time'

        state_variables = model.find_variable({'type': 'state'})
        assert len(state_variables) == 4

        import networkx as nx
        sorted_nodes = nx.lexicographical_topological_sort(graph, key=lambda x: str(x))

        sorted_nodes = list(sorted_nodes)
        assert str(sorted_nodes[0]) == '_environment$time'
        assert str(sorted_nodes[10]) == '_membrane$stim_period'
        assert str(sorted_nodes[20]) == '_sodium_channel$E_Na'
        assert str(sorted_nodes[-1]) == 'Derivative(_membrane$V, _environment$time)'

        for i, node in enumerate(sorted_nodes):
            print('%d. %r: %r' % (i, node, graph.nodes[node]['equation']))

        # use `dot -Tpng path.dot -o path.png`
        # nx.nx_agraph.write_dot(graph,
        #                        '/Users/tamuri/Desktop/path.dot')

        # free variable should not depend on anything
        time_dummy = free_variable[0]['sympy.Dummy']
        assert graph.in_degree(time_dummy) == 0

        # state variables should not depend on anything
        for variable in state_variables:
            dummy = variable['sympy.Dummy']
            assert graph.in_degree(dummy) == 0

        # check a node for dependencies
        dm_dt_node = sorted_nodes[29]
        assert str(dm_dt_node) == 'Derivative(_sodium_channel_m_gate$m, _environment$time)'
        assert 3 == graph.in_degree(dm_dt_node)
