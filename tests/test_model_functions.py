import os
from collections import OrderedDict

import pytest

import cellmlmanip
import cellmlmanip.rdf
from cellmlmanip.units import UnitStore

OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"


class TestModelFunctions():
    """ Tests for all functions on Model class """

    ###############################################################
    # fixtures

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

    @pytest.fixture(scope="class")
    def quantity_store(self, model):
        qs = UnitStore(model)
        for unit_name, unit_attributes in self.test_definitions.items():
            qs.add_custom_unit(unit_name, unit_attributes)
        return qs

    @pytest.fixture
    def model(scope='class'):
        return cellmlmanip.load_model(
            os.path.join(os.path.dirname(__file__), 'cellml_files', 'basic_ode.cellml'))

    @pytest.fixture
    def other_model(scope='class'):
        return cellmlmanip.load_model(
            os.path.join(os.path.dirname(__file__), 'cellml_files',
                         'aslanidi_model_2009.cellml'))

    ##########################################################
    # check equation graph property

    # also tested in test_hodgkin
    def test_graph_property(self, model):
        """ Tests that the graph property for Model has been constructed correctly. """

        graph1 = model.graph
        assert(len(graph1.nodes) == 3)
        names = ['env_ode$sv1', 'environment$time']
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

    #######################################################################
    # this section contains tests for each get_ function on Model

    # also tested in test_hodgkin
    def test_get_state_symbols(self, model):
        """ Tests Model.get_state_symbols() works correctly. """

        state_symbols = model.get_state_symbols()
        assert len(state_symbols) == 1
        assert state_symbols[0].name == 'env_ode$sv1'

    # also tested in test_hodgkin
    def test_get_free_variable_symbol(self, model):
        """ Tests Model.get_free_variable_symbol() works correctly. """

        free_variable_symbol = model.get_free_variable_symbol()
        assert free_variable_symbol.name == 'environment$time'

    def test_get_free_variable_symbol_1(self, other_model):
        """ Tests Model.get_free_variable_symbol() works correctly. """

        free_variable_symbol = other_model.get_free_variable_symbol()
        assert free_variable_symbol.name == 'environment$time'

    # also tested in test_aslanidi
    def test_get_initial_value(self, other_model):
        """ Tests Model.get_initial_value() works correctly. """

        membrane_voltage = other_model.get_symbol_by_ontology_term(OXMETA, "membrane_voltage")
        assert(other_model.get_initial_value(membrane_voltage) == -80.0)

    # also tested in test_hodgkin
    def test_get_derivative_symbols(self, model):
        """ Tests Model.get_derivative_symbols() works correctly. """

        derivs = model.get_derivative_symbols()
        assert len(derivs) == 1
        deriv = derivs[0]
        assert deriv.is_Derivative
        assert len(deriv.variables) == 1
        assert deriv.variables[0].is_Dummy
        assert deriv.variables[0].name == 'environment$time'

    # also tested by model_units
    def test_get_equations_for(self, model):
        """ Tests Model.get_equations_for() works correctly.
        Note: the basic model has no equations
        """

        symbol_a = model.get_symbol_by_cmeta_id("sv11")
        equation = model.get_equations_for([symbol_a])
        assert len(equation) == 0

    # also tested by model_units
    def test_get_equations_for_1(self, other_model):
        """ Tests Model.get_equations_for() works correctly. """

        symbol_a = other_model.get_symbol_by_ontology_term(OXMETA, "membrane_capacitance")
        equation = other_model.get_equations_for([symbol_a])
        assert len(equation) == 1
        assert equation[0].lhs == symbol_a
        assert equation[0].rhs == 5e-5

    def test_get_value(self, other_model):
        """ Tests Model.get_value() works correctly. """

        symbol_a = other_model.get_symbol_by_ontology_term(OXMETA, "membrane_capacitance")
        assert other_model.get_value(symbol_a) == 5e-5

    #################################################################
    # tests for get_symbol_ functions

    def test_get_symbol_by_cmeta_id(self, model):
        """ Tests Model.get_symbol_by_cmeta_id() works correctly. """

        sv11 = model.get_symbol_by_cmeta_id('sv11')
        assert sv11.name == 'env_ode$sv1'
        assert sv11.units == 'mV'

    def test_get_symbol_by_cmeta_id_2(self, other_model):
        """ Tests Model.get_symbol_by_cmeta_id() works correctly. """

        variable = other_model.get_symbol_by_cmeta_id('testcmeta')
        assert variable.name == 'intracellular_ion_concentrations$Na_i'
        assert variable.units == 'millimolar'

    def test_get_symbol_by_name(self, model):
        """ Tests Model.get_symbol_by_name() works correctly. """

        sv11 = model.get_symbol_by_name('env_ode$sv1')
        assert sv11.units == 'mV'

    def test_get_symbol_by_ontology_term(self, other_model):
        """ Tests Model.get_symbol_by_ontology_term() works correctly. """

        symbol_a = other_model.get_symbol_by_ontology_term(OXMETA, 'membrane_capacitance')
        assert symbol_a.name == 'membrane$Cm'
        assert symbol_a.units == 'nanoF'

    def test_get_symbols_by_rdf(self, other_model):
        """ Tests Model.get_symbols_by_rdf() works correctly. """

        symbol_a = other_model.get_symbols_by_rdf(('http://biomodels.net/biology-qualifiers/', 'is'),
                                                  (OXMETA, 'membrane_voltage'))
        assert len(symbol_a) == 1
        assert symbol_a[0].name == 'membrane$V'
        assert symbol_a[0].units == 'millivolt'

    ######################################################################
    # The functions listed for ontology/rdf are tested in test_rdf.py
    #
    # get_ontology_terms_by_symbol()
    # get_rdf_annotations()
    # get_rdf_value() - indirectly tested by test_get_rdf_annotations() in test_rdf.py
    # has_ontology_annotation()
