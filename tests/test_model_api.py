import os

import pytest
import sympy as sp

from cellmlmanip import parser
from cellmlmanip.model import NumberDummy


OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"


class TestModelAPI(object):
    """
    Tests various parts of the Model API.
    """

    @pytest.fixture(scope="class")
    def model(self):
        path = os.path.join(
            os.path.dirname(__file__), 'cellml_files', 'hodgkin_huxley_squid_axon_model_1952_modified.cellml')
        return parser.Parser(path).parse()

    @pytest.fixture(scope="class")
    def espinosa_model_1998_normal(self):
        path = os.path.join(
            os.path.dirname(__file__), 'cellml_files', 'espinosa_model_1998_normal.cellml')
        return parser.Parser(path).parse()

    def test_add_unit(self, model):
        """ Tests the add_unit method. """

        # TODO: Add unit tests for normal usage

        # Base units can't have attributes
        with pytest.raises(ValueError, match='can not be defined with unit attributes'):
            model.add_unit('unlikely_unit_name', [{'units': 'millivolt'}], base_units=True)

    def test_add_variable(self, model):
        """ Tests the add_variable method. """

        # TODO: Add unit tests for normal usage

        # Variable can't be added twice
        unit = 'millivolt'
        model.add_variable(name='varvar1', units=unit)
        model.add_variable(name='varvar2', units=unit)
        with pytest.raises(ValueError, match='already exists'):
            model.add_variable(name='varvar1', units=unit)

    def test_find_symbols_and_derivatives(self, model):
        """ Tests Model.find_symbols_and_derivatives. """

        # Test on single variable expressions
        t = model.get_free_variable_symbol()
        syms = model.find_symbols_and_derivatives([t])
        assert len(syms) == 1
        assert t in syms

        v = model.get_symbol_by_ontology_term(OXMETA, "membrane_voltage")
        dvdt = sp.Derivative(v, t)
        syms = model.find_symbols_and_derivatives([dvdt])
        assert len(syms) == 1
        assert sp.Derivative(v, t) in syms

        # Test on longer expressions
        x = sp.Float(1) + t * sp.sqrt(dvdt) - t
        syms = model.find_symbols_and_derivatives([x])
        assert len(syms) == 2
        assert t in syms
        assert sp.Derivative(v, t) in syms

        # Test on multiple expressions
        y = sp.Float(2) + v
        syms = model.find_symbols_and_derivatives([x, y])
        assert len(syms) == 3
        assert v in syms
        assert t in syms
        assert sp.Derivative(v, t) in syms

    def test_get_value(self, model):
        """ Tests Model.get_value() works correctly. """

        g = model.get_symbol_by_ontology_term(OXMETA, 'membrane_fast_sodium_current_conductance')
        assert model.get_value(g) == 120

    def test_graph_for_dae(self):
        """ Checks if writing a DAE in a model raises an exceptions. """

        # Parsing should be OK
        path = os.path.join(os.path.dirname(__file__), 'cellml_files', '4.algebraic_ode_model.cellml')
        model = parser.Parser(path).parse()

        # But equation graph will raise error (if accessed)
        with pytest.raises(RuntimeError, match='DAEs are not supported'):
            model.graph

    def test_set_equation(self, model):
        """ Tests replacing an equation in a model. """

        # Get model, assert that V is a state variable
        v = model.get_symbol_by_ontology_term(OXMETA, 'membrane_voltage')
        assert v.type == 'state'

        # Now clamp it to -80mV
        rhs = model.add_number(-80, str(v.units))
        model.set_equation(v, rhs)

        # Check that V is no longer a state
        v = model.get_symbol_by_ontology_term(OXMETA, 'membrane_voltage')
        assert v.type != 'state'

        # TODO: Get dvdt_unit in a more sensible way
        # See: https://github.com/ModellingWebLab/cellmlmanip/issues/133

        # Now make V a state again
        t = model.get_symbol_by_ontology_term(OXMETA, 'time')
        lhs = sp.Derivative(v, t)
        dvdt_units = 'unlikely_unit_name'
        model.add_unit(dvdt_units, [
            {'units': str(v.units)},
            {'units': str(t.units), 'exponent': -1},
        ])
        rhs = model.add_number(0, dvdt_units)
        model.set_equation(lhs, rhs)

        # Check that V is a state again
        v = model.get_symbol_by_ontology_term(OXMETA, 'membrane_voltage')
        assert v.type == 'state'

        # Set equation for a newly created variable
        lhs = model.add_variable(name='an_incredibly_unlikely_variable_name', units=str(v.units))
        rhs = model.add_number(12, str(v.units))
        model.set_equation(lhs, rhs)

    def test_get_equations_for(self, espinosa_model_1998_normal):
        model = espinosa_model_1998_normal

        eqs_no_recurse = model.get_equations_for(model.get_derivative_symbols(), recurse=False)
        eqs_no_recurse_strip_units = \
            model.get_equations_for(model.get_derivative_symbols(), recurse=False, strip_units=False)
        assert len(eqs_no_recurse) == len(eqs_no_recurse_strip_units) == 80
        assert isinstance(eqs_no_recurse[0].rhs, sp.numbers.Float)
        assert isinstance(eqs_no_recurse_strip_units[0].rhs, NumberDummy)
        print(type(eqs_no_recurse_strip_units[0].rhs))

        eqs = model.get_equations_for(model.get_derivative_symbols())
        eqs2 = model.get_equations_for(model.get_derivative_symbols(), keep_unused_eqs=False)
        assert (len(eqs)) == (len(eqs2)) + 1 == 146
        diff = [eq for eq in eqs if eq not in eqs2]
        assert len(diff) == 1
        assert str(diff[0].lhs) == 'calcium_release$VoltDep'
