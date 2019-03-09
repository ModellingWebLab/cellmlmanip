import os

import networkx as nx
import pytest
import sympy

from cellmlmanip import load_model


class TestHodgkin:
    @pytest.fixture(scope="class")
    def model(self):
        hodgkin_cellml = os.path.join(
            os.path.dirname(__file__),
            "cellml_files",
            "hodgkin_huxley_squid_axon_model_1952_modified.cellml"
        )
        return load_model(hodgkin_cellml)

    @pytest.fixture(scope="class")
    def graph(self, model):
        return model.get_equation_graph()

    def test_counts(self, model):
        # https://models.cellml.org/exposure/5d116522c3b43ccaeb87a1ed10139016/hodgkin_huxley_1952_variant01.cellml/cellml_math
        assert len(model.equations) == 17

    def test_connections(self, model):
        target = model.get_dummy_data('sodium_channel$h')
        source = model.get_dummy_data('sodium_channel_h_gate$h')
        assert target.assigned_to == source.dummy

    def test_equation_units(self, model):
        equation = model.equations[2]
        lhs_units = model.units.summarise_units(equation.lhs)
        assert lhs_units == model.units.ureg.millivolt

    def test_check_left_right_units(self, model):
        for e in model.equations:
            model.check_left_right_units_equal(e)

    def test_equation_graph(self, graph, model):
        assert len(graph.nodes) == 32

        free_variable = model.find_variable({'type': 'free'})
        assert len(free_variable) == 1
        free_variable = free_variable[0]
        assert free_variable.cmeta_id == 'time'
        assert graph.node[free_variable.dummy]['variable_type'] == 'free'
        assert free_variable.cmeta_id == graph.node[free_variable.dummy]['cmeta_id']

        state_variables = model.find_variable({'type': 'state'})
        assert len(state_variables) == 4
        state_variable = state_variables[0]
        assert graph.node[state_variable.dummy]['variable_type'] == 'state'
        assert state_variable.cmeta_id == graph.node[state_variable.dummy]['cmeta_id']

        sorted_nodes = nx.lexicographical_topological_sort(graph, key=str)

        sorted_nodes = list(sorted_nodes)
        assert sorted_nodes[0].name == 'environment$time'
        assert sorted_nodes[10].name == 'membrane$stim_period'
        assert sorted_nodes[20].name == 'sodium_channel$E_Na'
        assert str(sorted_nodes[-1]) == 'Derivative(_membrane$V, _environment$time)'

        # check all cmeta ids have been added
        for node in sorted_nodes:
            # derivative nodes depend on state and free variable nodes
            if not node.is_Derivative:
                variable = model.get_dummy_data(node)
                for key in ['cmeta_id', 'name']:
                    if getattr(variable, key, None):
                        assert getattr(variable, key) == graph.nodes[node][key]

                # only state variables should have initial_values
                if graph.nodes[node].get('variable_type', '') == 'state':
                    assert (float(variable.initial_value) ==
                            float(graph.nodes[node]['initial_value']))
                else:
                    assert 'initial_value' not in graph.nodes[node]

        # for i, node in enumerate(sorted_nodes):
        #     print('%d. %r: %r' % (i, node, graph.nodes[node]['equation']))

        # use `dot -Tpng path.dot -o path.png`
        # nx.nx_agraph.write_dot(graph,
        #                        '/Users/tamuri/Desktop/path.dot')

        # free variable should not depend on anything
        time_dummy = free_variable.dummy
        assert graph.in_degree(time_dummy) == 0

        # state variables should not depend on anything
        for variable in state_variables:
            dummy = variable.dummy
            assert graph.in_degree(dummy) == 0

        # check a node for dependencies
        dm_dt_node = sorted_nodes[29]
        assert str(dm_dt_node) == 'Derivative(_sodium_channel_m_gate$m, _environment$time)'
        assert 3 == graph.in_degree(dm_dt_node)

        # check attributes on membrane capacitance
        membrane_Cm = sorted_nodes[2]
        assert membrane_Cm.name == 'membrane$Cm'
        assert graph.node[membrane_Cm]['cmeta_id'] == 'membrane_capacitance'
        assert graph.node[membrane_Cm]['variable_type'] == 'parameter'

    def test_derivative_symbols(self, model):
        """
        Derivative(_membrane$V, _environment$time)
        Derivative(_sodium_channel_m_gate$m, _environment$time)
        Derivative(_sodium_channel_h_gate$h, _environment$time)
        Derivative(_potassium_channel_n_gate$n, _environment$time)
        """
        derivatives = model.get_derivative_symbols()
        assert len(derivatives) == 4

    def test_state_symbols(self, model):
        """
        _membrane$V
        _sodium_channel_m_gate$m
        _sodium_channel_h_gate$h
        _potassium_channel_n_gate$n
        """
        state_symbols = model.get_state_symbols()
        assert len(state_symbols) == 4

    def test_free_variable_symbol(self, model):
        """
        _environment$time
        """
        free_variable_symbol = model.get_free_variable_symbol()
        assert free_variable_symbol.name == 'environment$time'

    def test_evaluation(self, graph, model):
        """
        From Michael:

        These are my calculations for the derivatives of the state variables
        in the hodgkin-huxley model, evaluated from the initial state (so the
        initial values of the state variables)

        Reading model from hodgkin_huxley_squid_axon_model_1952_modified.mmt...
        Model hodgkin_huxley_squid_axon_model_1952_modified read successfully.
        Evaluating state vector derivatives...

        -------------------------------------------------------------------------------
        Name                        Initial value             Derivative at t=0
        -------------------------------------------------------------------------------
        membrane.V                  -7.50000000000000000e+01  -6.00768750000000740e-01
        sodium_channel_m_gate.m      5.00000000000000028e-02   1.23855383553985177e-02
        sodium_channel_h_gate.h      5.99999999999999978e-01  -4.55523906540064583e-04
        potassium_channel_n_gate.n   3.25000000000000011e-01  -1.34157228632045961e-03
        -------------------------------------------------------------------------------```

        The free variable time is set to 0

        Numbers are converted to string with `'{:<1.17e}'.format(number)`,
        which I think should give the maximum precision for any double precision float

        (as in, converting to and back from string shouldn't lose
        accuracy, it does _not_ mean all printed digits are accurate)
        """

        # initial values for the free and state variables
        initials = {
            'environment$time': 0.0,
            'membrane$V': -7.50000000000000000e+01,
            'sodium_channel_m_gate$m': 5.00000000000000028e-02,
            'sodium_channel_h_gate$h': 5.99999999999999978e-01,
            'potassium_channel_n_gate$n': 3.25000000000000011e-01
        }

        # the calculated derivatives at the initial conditions
        evaluated_derivatives = {
            'membrane$V': -6.00768750000000740e-01,
            'sodium_channel_m_gate$m': 1.23855383553985177e-02,
            'sodium_channel_h_gate$h': -4.55523906540064583e-04,
            'potassium_channel_n_gate$n': -1.34157228632045961e-03
        }

        sorted_symbols = nx.lexicographical_topological_sort(graph, key=str)

        # collects all the evaluated lh-sides of the equations
        evaluated_symbols = dict()

        # saves the evaluated lh-side for the state variables
        state_derivatives = dict()

        # loop over each of the equations in the model (topologically sorted)
        for symbol in sorted_symbols:
            equation = graph.nodes[symbol]['equation']

            # if this symbol is calculated using an equation
            if equation is not None:

                # substitute all symbols in the rhs of the equation
                eq_substituted = equation.rhs.subs(evaluated_symbols)

                # if any symbols remain, we have a problem
                remaining_symbols = eq_substituted.atoms(sympy.Symbol)
                if remaining_symbols:
                    for remaining_symbol in remaining_symbols:
                        pytest.fail("Unresolved symbol %s in %s" % (remaining_symbol, equation))

                # calculate the result
                eq_evaluated = eq_substituted.evalf()

                # save the calculation for this symbol (to be used in the next equations)
                evaluated_symbols[equation.lhs] = eq_evaluated

                # if the symbol on the lhs is a derivative
                if equation.lhs.is_Derivative:
                    # save the calculation for testing
                    state_derivatives[equation.lhs.free_symbols.pop()] = eq_evaluated

            # otherwise the symbol doesn't have an equation
            # if the symbol has an initial value
            elif symbol.name in initials:
                # add the symbol's initial value to the substitution dictionary
                evaluated_symbols[symbol] = sympy.Number(initials[symbol.name])
            else:
                # something has gone wrong - no equation or initial value
                pytest.fail("Symbol " + str(symbol) + " does not have equation or initial value.")

        # for each calculated state variable
        for state_symbol, evaluated_deriv in state_derivatives.items():
            # check evaluation against expected
            expected = evaluated_derivatives[state_symbol.name]
            actual = evaluated_deriv
            assert float(actual) == pytest.approx(expected)
