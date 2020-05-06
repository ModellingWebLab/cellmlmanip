import networkx
import pytest
import sympy


def test_model_evaluation(hh_model):
    """
    Tests the calculated derivatives of a single model match a known reference.
    """
    graph = hh_model.graph_with_sympy_numbers

    # Initial values for the free and state variables
    initial_values = {
        'environment$time': 0.0,
        'membrane$V': -75,
        'sodium_channel_m_gate$m': 0.05,
        'sodium_channel_h_gate$h': 0.6,
        'potassium_channel_n_gate$n': 0.325,
    }

    # Evaluated derivatives at the initial conditions
    evaluated_derivatives = {
        'membrane$V': -6.00768750000000740e-01,
        'sodium_channel_m_gate$m': 1.23855383553985177e-02,
        'sodium_channel_h_gate$h': -4.55523906540064583e-04,
        'potassium_channel_n_gate$n': -1.34157228632045961e-03,
    }

    # collects all the evaluated lh-sides of the equations
    evaluated_symbols = dict()

    # saves the evaluated lh-side for the state variables
    state_derivatives = dict()

    # loop over each of the equations in the model (topologically sorted)
    sorted_symbols = networkx.lexicographical_topological_sort(graph, key=str)
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
        elif symbol.name in initial_values:
            # add the symbol's initial value to the substitution dictionary
            evaluated_symbols[symbol] = sympy.Number(initial_values[symbol.name])
        else:
            # something has gone wrong - no equation or initial value
            pytest.fail("Symbol " + str(symbol) + " does not have equation or initial value.")

    # Check evaluations against expected values
    for state_symbol, evaluated_deriv in state_derivatives.items():
        expected = evaluated_derivatives[state_symbol.name]
        actual = evaluated_deriv
        assert float(actual) == pytest.approx(expected)
