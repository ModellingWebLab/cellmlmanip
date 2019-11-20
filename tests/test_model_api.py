import os
import pytest
import sympy as sp

from cellmlmanip import parser


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

    def test_dae(self):
        """ Checks if writing a DAE in a model raises an exceptions. """

        # Parsing should be OK
        path = os.path.join(os.path.dirname(__file__), 'cellml_files', '4.algebraic_ode_model.cellml')
        model = parser.Parser(path).parse()

        # But equation graph will raise error
        with pytest.raises(RuntimeError, match='DAEs are not supported'):
            model.get_equation_graph()

    def test_replace_equation(self, model):
        """ Tests replacing an equation in a model. """

        # Get model, assert that V is a state variable
        model.get_equation_graph()
        v = model.get_symbol_by_ontology_term(OXMETA, 'membrane_voltage')
        v_meta = model.get_meta_dummy(v)
        assert v_meta.type == 'state'

        # Now clamp it to -80mV
        v_unit = v_meta.units
        rhs = model.add_number(number=sp.Number(-80), units=str(v_unit))
        model.replace_equation(sp.Eq(v, rhs))

        # Check that V is no longer a state
        model.get_equation_graph()
        v = model.get_symbol_by_ontology_term(OXMETA, 'membrane_voltage')
        v_meta = model.get_meta_dummy(v)
        assert v_meta.type != 'state'

        # TODO: Get dvdt_unit in a more sensible way
        # See: https://github.com/ModellingWebLab/cellmlmanip/issues/133

        # Now make V a state again
        t = model.get_symbol_by_ontology_term(OXMETA, 'time')
        t_meta = model.get_meta_dummy(t)
        t_unit = t_meta.units
        lhs = sp.Derivative(v, t)
        dvdt_unit = 'unlikely_unit_name'
        model.add_unit(dvdt_unit, [
            {'units': str(v_unit)},
            {'units': str(t_unit), 'exponent': -1},
        ])
        rhs = model.add_number(number=sp.Number(0), units=dvdt_unit)
        model.replace_equation(sp.Eq(lhs, rhs))

        # Check that V is a state again
        model.get_equation_graph()
        v = model.get_symbol_by_ontology_term(OXMETA, 'membrane_voltage')
        v_meta = model.get_meta_dummy(v)
        assert v_meta.type == 'state'
