import logging
import os

import networkx as nx
import pytest
import sympy
from sympy.physics.units import Quantity

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

    @pytest.fixture(scope="class")
    def graph(self, model):
        return model.get_equation_graph()

    def test_counts(self, model):
        # https://models.cellml.org/exposure/5d116522c3b43ccaeb87a1ed10139016/hodgkin_huxley_1952_variant01.cellml/cellml_math
        assert len(model.components) == 8
        eq_count = len(list(model.equations))
        assert eq_count == 17

    def test_setup_connections(self, model):
        model.make_connections()
        target = model.components['sodium_channel'].variables['h']
        source = model.components['sodium_channel_h_gate'].variables['h']
        assert target['assignment'] == source['sympy.Dummy']

    def test_add_units_to_equations(self, model):
        model.add_units_to_equations()
        equation = model.components['sodium_channel'].equations[0]
        lhs_units = model.units.summarise_units(equation.lhs)
        assert lhs_units == model.units.store['millivolt']

    def test_check_left_right_units(self, model):
        for e in model.equations:
            model.check_left_right_units_equal(e)

    def test_get_equations(self, graph, model):
        assert len(graph.nodes) == 32

        free_variable = model.find_variable({'type': 'free'})
        assert len(free_variable) == 1
        free_variable = free_variable[0]
        assert free_variable['cmeta:id'] == 'time'
        assert graph.node[free_variable['sympy.Dummy']]['variable_type'] == 'free'
        assert free_variable['cmeta:id'] == graph.node[free_variable['sympy.Dummy']]['cmeta:id']

        state_variables = model.find_variable({'type': 'state'})
        assert len(state_variables) == 4
        state_variable = state_variables[0]
        assert graph.node[state_variable['sympy.Dummy']]['variable_type'] == 'state'
        assert state_variable['cmeta:id'] == graph.node[state_variable['sympy.Dummy']]['cmeta:id']

        sorted_nodes = nx.lexicographical_topological_sort(graph, key=lambda x: str(x))

        sorted_nodes = list(sorted_nodes)
        assert str(sorted_nodes[0]) == '_environment$time'
        assert str(sorted_nodes[10]) == '_membrane$stim_period'
        assert str(sorted_nodes[20]) == '_sodium_channel$E_Na'
        assert str(sorted_nodes[-1]) == 'Derivative(_membrane$V, _environment$time)'

        # check all cmeta ids have been added
        for node in sorted_nodes:
            # derivative nodes depend on state and free variable nodes
            if not node.is_Derivative:
                variable = model.find_variable({'sympy.Dummy': node})
                assert len(variable) == 1
                for key in ['cmeta:id', 'name']:
                    if key in variable[0]:
                        assert variable[0][key] == graph.nodes[node][key]
                if 'initial_value' in variable[0]:
                    assert float(variable[0]['initial_value']) == float(graph.nodes[node]['initial_value'])

        # for i, node in enumerate(sorted_nodes):
        #     print('%d. %r: %r' % (i, node, graph.nodes[node]['equation']))

        # use `dot -Tpng path.dot -o path.png`
        # nx.nx_agraph.write_dot(graph,
        #                        '/Users/tamuri/Desktop/path.dot')

        # free variable should not depend on anything
        time_dummy = free_variable['sympy.Dummy']
        assert graph.in_degree(time_dummy) == 0

        # state variables should not depend on anything
        for variable in state_variables:
            dummy = variable['sympy.Dummy']
            assert graph.in_degree(dummy) == 0

        # check a node for dependencies
        dm_dt_node = sorted_nodes[29]
        assert str(dm_dt_node) == 'Derivative(_sodium_channel_m_gate$m, _environment$time)'
        assert 3 == graph.in_degree(dm_dt_node)

        # check that units have been added to parameters
        leakage_var = sorted_nodes[1]
        assert str(leakage_var) == '_leakage_current$g_L'
        lcv_equation = graph.node[leakage_var]['equation']
        rhs_unit = model.units.summarise_units(lcv_equation.rhs)
        lhs_unit = model.units.summarise_units(lcv_equation.lhs)
        assert rhs_unit == model.units.store['milliS_per_cm2']
        assert rhs_unit == lhs_unit

        # check attributes on membrane capacitance
        membrane_Cm = sorted_nodes[2]
        assert str(membrane_Cm) == '_membrane$Cm'
        assert graph.node[membrane_Cm]['cmeta:id'] == 'membrane_capacitance'
        assert graph.node[membrane_Cm]['variable_type'] == 'parameter'

    def test_evaluation(self, graph):
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
            '_environment$time': 0.0,
            '_membrane$V': -7.50000000000000000e+01,
            '_sodium_channel_m_gate$m': 5.00000000000000028e-02,
            '_sodium_channel_h_gate$h': 5.99999999999999978e-01,
            '_potassium_channel_n_gate$n': 3.25000000000000011e-01
        }

        # the calculated derivatives at the initial conditions
        evaluated_derivatives = {
            '_membrane$V': -6.00768750000000740e-01,
            '_sodium_channel_m_gate$m': 1.23855383553985177e-02,
            '_sodium_channel_h_gate$h': -4.55523906540064583e-04,
            '_potassium_channel_n_gate$n': -1.34157228632045961e-03
        }

        sorted_symbols = nx.lexicographical_topological_sort(graph, key=lambda x: str(x))

        def __remove_quantities(eq_with_quantities):
            """Replaces all quantity symbols in the equation with 1"""
            quantities = {q: 1 for q in eq_with_quantities.atoms(Quantity)}
            eq_without_quantities = eq_with_quantities.subs(quantities, simultaneous=True)
            return eq_without_quantities

        # collects all the evaluated lh-sides of the equations
        evaluated_symbols = dict()

        # saves the evaluated lh-side for the state variables
        state_derivatives = dict()

        # loop over each of the equations in the model (topologically sorted)
        for symbol in sorted_symbols:
            equation = graph.nodes[symbol]['equation']

            # if this symbol is calculated using an equation
            if equation is not None:

                # remove all quantity symbols from the equation
                eq_no_units = __remove_quantities(equation.rhs)

                # substitute all symbols in the rhs of the equation
                eq_substituted = eq_no_units.subs(evaluated_symbols)

                # if any symbols remain, we have a problem
                remaining_symbols = eq_substituted.atoms(sympy.Symbol)
                if remaining_symbols:
                    for remaining_symbol in remaining_symbols:
                        # TODO: we need to handle 0.0 - best place to put this?
                        if str(remaining_symbol) == '_0.0':
                            eq_substituted = eq_substituted.subs({remaining_symbol: 0.0})
                        else:
                            pytest.fail("Unresolved symbol" + remaining_symbol + " in " + equation)

                # calculate the result
                eq_evaluated = eq_substituted.evalf()

                # save the calculation for this symbol (to be used in the next equations)
                lhs_without_units = __remove_quantities(equation.lhs)
                evaluated_symbols[lhs_without_units] = eq_evaluated

                # if the symbol on the lhs is a derivative
                if lhs_without_units.is_Derivative:
                    # save the calculation for testing
                    state_derivatives[lhs_without_units.free_symbols.pop()] = eq_evaluated

            # otherwise the symbol doesn't have an equation
            # if the symbol has an initial value
            elif str(symbol) in initials:
                # add the symbol's initial value to the substitution dictionary
                evaluated_symbols[symbol] = sympy.Number(initials[str(symbol)])
            else:
                # something has gone wrong - no equation or initial value
                pytest.fail("Symbol " + str(symbol) + " does not have equation or initial value.")

        # for each calculated state variable
        for state_symbol, evaluated_deriv in state_derivatives.items():
            # check evaluation against expected
            expected = evaluated_derivatives[str(state_symbol)]
            actual = evaluated_deriv
            assert float(actual) == pytest.approx(expected)
