import os
from collections import OrderedDict

import pytest
import sympy as sp
import cellmlmanip
import cellmlmanip.rdf

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
    # this section contains tests for each get_XXX function on Model

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
        Note: the basic model has 1 equation dsv11/dt = 1 which is not
        related to a symbol and so has no equations that can be retrieved
        by symbol
        """

        symbol_a = model.get_symbol_by_cmeta_id("sv11")
        equation = model.get_equations_for([symbol_a])
        assert len(equation) == 0

        symbol_t = model.get_symbol_by_name("environment$time")

        equation1 = model.get_equations_for([symbol_t])
        assert len(equation1) == 0

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
    # tests for get_symbol_XXX functions

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

    #########################################################################
    # test add functions
    # note some tests for fail cases are repeated in TestModelAPI class

    def test_add_equation(self, model):
        """ Tests the Model.add_equation method.
        """
        assert len(model.equations) == 1
        # so we are adding
        # newvar2 = newvar + newvar1
        # but need to also add newvar1 = 2; newvar = 2 in order or the graph to resolve correctly
        model.add_variable(name='newvar', units='mV')
        symbol = model.get_symbol_by_name('newvar')
        model.add_variable(name='newvar1', units='mV')
        symbol1 = model.get_symbol_by_name('newvar1')
        model.add_variable(name='newvar2', units='mV')
        symbol2 = model.get_symbol_by_name('newvar2')
        model.add_equation(sp.Eq(symbol, 2.0))
        model.add_equation(sp.Eq(symbol1, 2.0))
        model.add_equation(sp.Eq(symbol2, sp.Add(symbol, symbol1)))
        assert len(model.equations) == 4
        eqn = model.get_equations_for([symbol2])
        assert len(eqn) == 3
        assert eqn[0].lhs == symbol
        assert eqn[0].rhs == 2.0
        assert eqn[1].lhs == symbol1
        assert eqn[1].rhs == 2.0
        assert eqn[2].lhs == symbol2
        assert eqn[2].rhs == sp.Add(symbol, symbol1)

    def test_set_equation(self, model):
        """ Tests the Model.set_equation method.
        """
        assert len(model.equations) == 1
        # so we are adding
        # newvar2 = newvar + newvar1
        # but need to also add newvar1 = 2; newvar = 2 in order or the graph to resolve correctly
        model.add_variable(name='newvar', units='mV')
        symbol = model.get_symbol_by_name('newvar')
        model.add_variable(name='newvar1', units='mV')
        symbol1 = model.get_symbol_by_name('newvar1')
        model.add_variable(name='newvar2', units='mV')
        symbol2 = model.get_symbol_by_name('newvar2')
        model.set_equation(symbol, 2.0)
        model.set_equation(symbol1, 2.0)
        model.set_equation(symbol2, sp.Add(symbol, symbol1))
        assert len(model.equations) == 4
        eqn = model.get_equations_for([symbol2])
        assert len(eqn) == 3
        assert eqn[0].lhs == symbol
        assert eqn[0].rhs == 2.0
        assert eqn[1].lhs == symbol1
        assert eqn[1].rhs == 2.0
        assert eqn[2].lhs == symbol2
        assert eqn[2].rhs == sp.Add(symbol, symbol1)

    def test_add_number(self, model):
        """ Tests the Model.add_number method. """
        number2 = model.add_number(2.0, 'mV')
        assert number2.is_Dummy

    def test_add_unit(self, model):
        """ Tests the Model.add_unit method. """
        assert len(model.units.custom_defined) == 5

        assert 'uF' not in model.units.custom_defined
        model.add_unit('uF', [{'units': 'farad', 'prefix': 'micro'}])
        assert len(model.units.custom_defined) == 6
        assert 'uF' in model.units.custom_defined

        # repeated in TestModelAPI
        # Base units can't have attributes
        with pytest.raises(ValueError, match='can not be defined with unit attributes'):
            model.add_unit('unlikely_unit_name', [{'units': 'millivolt'}], base_units=True)
        assert len(model.units.custom_defined) == 6

    def test_add_variable(self, model):
        """ Tests the Model.add_variable() method. """

        assert len(model.variables()) == 3
        with pytest.raises(KeyError):
            model.get_symbol_by_name('newvar') is None

        model.add_variable(name='newvar', units='mV')
        assert len(model.variables()) == 4
        assert model.get_symbol_by_name('newvar')

        # Repeatd in TestModelAPI
        # Variable can't be added twice
        unit = 'mV'
        model.add_variable(name='varvar1', units=unit)
        model.add_variable(name='varvar2', units=unit)
        assert len(model.variables()) == 6
        with pytest.raises(ValueError, match='already exists'):
            model.add_variable(name='varvar1', units=unit)

    # TO DO
    def test_add_rdf(self, model):
        pass

    ###################################################################
    # this section is for other functions

    # TO DO
    def test_connect_variables(selfself, model):
        pass

