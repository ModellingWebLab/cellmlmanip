import os
from collections import OrderedDict

import pytest

import cellmlmanip
import cellmlmanip.rdf
from cellmlmanip.units import UnitStore

OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"


class TestModelFunctions():
    # These represent CellML <units><unit>...</unit></units> elements
    test_definitions = OrderedDict({
        'ms': [{'units': 'second', 'prefix': 'milli'}],
        'usec': [{'units': 'second', 'prefix': 'micro'}],
        'mV': [{'units': 'volt', 'prefix': 'milli'}],
        'uV': [{'units': 'volt', 'prefix': 'micro'}],
        'mM': [{'prefix': 'milli', 'units': 'mole'}, {'units': 'litre', 'exponent': '-1'}],
        'milli_mole': [{'prefix': 'milli', 'units': 'mole'}],
        'millisecond': [{'prefix': 'milli', 'units': 'second'}],
    })
    test_definitions.update({
        'per_ms': [{'units': 'ms', 'exponent': '-1'}],
        'per_mV': [{'units': 'volt', 'prefix': 'milli', 'exponent': '-1'}],
        'mV_per_ms': [{'units': 'mV', 'exponent': '1'}, {'units': 'ms', 'exponent': '-1'}],
        'mV_per_s': [{'units': 'mV', 'exponent': '1'}, {'units': 'second', 'exponent': '-1'}],
        'mV_per_usec': [
            {'units': 'mV', 'exponent': '1'}, {'prefix': 'micro', 'units': 'second', 'exponent': '-1'}],
        'mM_per_ms': [{'units': 'mM'}, {'units': 'ms', 'exponent': '-1'}],
        'ms_power_prefix': [{'prefix': '-3', 'units': 'second'}],
        'ms_with_multiplier': [{'multiplier': 0.001, 'units': 'second'}],
    })

    @pytest.fixture
    def model(scope='class'):
        return cellmlmanip.load_model(
            os.path.join(os.path.dirname(__file__), 'cellml_files', "beeler_reuter_model_1977.cellml"))

    @pytest.fixture
    def basic_model(scope='class'):
        return cellmlmanip.load_model(
            os.path.join(os.path.dirname(__file__), 'cellml_files', "basic_ode.cellml"))

    @pytest.fixture(scope="class")
    def graph(self, model):
        return model.get_equation_graph()

    @pytest.fixture(scope="class")
    def quantity_store(self, model):
        qs = UnitStore(model)
        for unit_name, unit_attributes in self.test_definitions.items():
            qs.add_custom_unit(unit_name, unit_attributes)
        return qs

    # also tested in test_hodgkin
    def test_get_state_symbols(self, model):
        state_symbols = model.get_state_symbols()
        assert len(state_symbols) == 8

    # also tested in test_hodgkin
    def test_get_derivative_symbols(self, model):
        derivs = model.get_derivative_symbols()
        assert len(derivs) == 8

    # also tested in test_hodgkin
    def test_get_free_variable_symbol(self, model):
        free_variable_symbol = model.get_free_variable_symbol()
        assert free_variable_symbol.name == 'environment$time'

    # also tested in test_aslanidi
    def test_get_initial_value(self, model):
        membrane_voltage = model.get_symbol_by_ontology_term(OXMETA, "membrane_voltage")
        assert(model.get_initial_value(membrane_voltage) == -84.624)

    # also tested in test_hodgkin
    def test_get_equation_graph(self, basic_model):
        graph1 = basic_model.get_equation_graph(True)
        assert(len(graph1.nodes) == 2)
        names = ['ode$sv1', 'ode$time']
        for v in graph1:
            if not v.is_Derivative:
                assert (v.name in names)
            else:
                for a in v._args:
                    if a.is_Dummy:
                        assert(a.name in names)
                    else:
                        for b in a._args:
                            if b.is_Dummy:
                                assert (b.name in names)

    # also tested by model_units
    def test_get_equations_for(self, basic_model):
        basic_model.get_equation_graph(True)
        symbol_a = basic_model.get_symbol_by_cmeta_id("sv11")
        equation = basic_model.get_equations_for([symbol_a])
        assert len(equation) == 1
        assert equation[0].lhs == symbol_a
        assert equation[0].rhs == 2.0
