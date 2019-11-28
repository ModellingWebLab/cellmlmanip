import os
from collections import OrderedDict

import pytest
import sympy as sp

from cellmlmanip import parser
from cellmlmanip.model import VariableDummy

from . import shared


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
    def local_model(scope='function'):
        """ Fixture to load a local copy of  the basic_ode model that may get modified. """
        return shared.load_model('basic_ode')

    @pytest.fixture
    def local_hh_model(scope='function'):
        """ Fixture to load a local copy of  the hodgkin_huxley_squid_axon_model_1952_modified
        model that may get modified. """
        return shared.load_model('hodgkin_huxley_squid_axon_model_1952_modified')

    ##########################################################
    # check equation graph property

    # also tested in test_hodgkin
    def test_graph_property(self, basic_model):
        """ Tests that the graph property for Model has been constructed correctly. """

        graph1 = basic_model.graph
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

    def test_graph_for_dae(self):
        """ Checks if writing a DAE in a model raises an exceptions. """

        # Parsing should be OK
        path = os.path.join(os.path.dirname(__file__), 'cellml_files', '4.algebraic_ode_model.cellml')
        model = parser.Parser(path).parse()

        # But equation graph will raise error (if accessed)
        with pytest.raises(RuntimeError, match='DAEs are not supported'):
            model.graph

    #######################################################################
    # this section contains tests for each get_XXX function on Model

    def test_get_state_symbols(self, basic_model):
        """ Tests Model.get_state_symbols() works correctly. """

        state_symbols = basic_model.get_state_symbols()
        assert len(state_symbols) == 1
        assert state_symbols[0].name == 'env_ode$sv1'

    def test_get_state_symbols2(self, aslanidi_model):
        """ Tests Model.get_state_symbols() works correctly. """
        state_symbols = aslanidi_model.get_state_symbols()
        assert len(state_symbols) == 29
        assert str(state_symbols) == \
            '[_membrane$V, _sodium_current_m_gate$m, _sodium_current_h1_gate$h1, _sodium_current_h2_gate$h2, '\
            '_L_type_Ca_channel_d_L_gate$d_L, _L_type_Ca_channel_f_L_gate$f_L, '\
            '_T_type_Ca_channel_d_T_gate$d_T, _T_type_Ca_channel_f_T_gate$f_T, '\
            '_Ca_independent_transient_outward_K_current_r_gate$r, '\
            '_Ca_independent_transient_outward_K_current_s1_gate$s1, '\
            '_Ca_independent_transient_outward_K_current_s2_gate$s2, '\
            '_Ca_independent_transient_outward_K_current_s3_gate$s3, _delayed_rectifier_K_current_z_gate$z, '\
            '_delayed_rectifier_K_current_pa_gate$p_a, _delayed_rectifier_K_current_pi_gate$p_i, '\
            '_intracellular_ion_concentrations$Na_i, _intracellular_ion_concentrations$Ca_i, '\
            '_intracellular_ion_concentrations$K_i, _intracellular_Ca_buffering$O_C, '\
            '_intracellular_Ca_buffering$O_TC, _intracellular_Ca_buffering$O_TMgC, '\
            '_intracellular_Ca_buffering$O_TMgMg, _cleft_space_ion_concentrations$K_c, '\
            '_Ca_handling_by_the_SR$Ca_rel, _Ca_handling_by_the_SR$Ca_up, _Ca_handling_by_the_SR$O_Calse, '\
            '_Ca_handling_by_the_SR$F1, _Ca_handling_by_the_SR$F2, _Ca_handling_by_the_SR$F3]'

    # also tested in test_hodgkin
    def test_get_free_variable_symbol(self, basic_model):
        """ Tests Model.get_free_variable_symbol() works correctly. """

        free_variable_symbol = basic_model.get_free_variable_symbol()
        assert free_variable_symbol.name == 'environment$time'

    def test_get_free_variable_symbol_1(self, aslanidi_model):
        """ Tests Model.get_free_variable_symbol() works correctly. """

        free_variable_symbol = aslanidi_model.get_free_variable_symbol()
        assert free_variable_symbol.name == 'environment$time'

    # also tested in test_aslanidi
    def test_get_initial_value(self, aslanidi_model, OXMETA):
        """ Tests Model.get_initial_value() works correctly. """

        membrane_voltage = aslanidi_model.get_symbol_by_ontology_term(OXMETA, "membrane_voltage")
        assert(aslanidi_model.get_initial_value(membrane_voltage) == -80.0)

    def test_get_derivative_symbols(self, basic_model):
        """ Tests Model.get_derivative_symbols() works correctly. """

        derivs = basic_model.get_derivative_symbols()
        assert len(derivs) == 1
        deriv = derivs[0]
        assert deriv.is_Derivative
        assert len(deriv.variables) == 1
        assert deriv.variables[0].is_Dummy
        assert deriv.variables[0].name == 'environment$time'

    def test_get_derivative_symbols2(self, aslanidi_model):
        derivs = aslanidi_model.get_derivative_symbols()
        assert len(derivs) == 29

        assert str(derivs) == '[Derivative(_membrane$V, _environment$time), '\
            'Derivative(_sodium_current_m_gate$m, _environment$time), '\
            'Derivative(_sodium_current_h1_gate$h1, _environment$time), '\
            'Derivative(_sodium_current_h2_gate$h2, _environment$time), '\
            'Derivative(_L_type_Ca_channel_d_L_gate$d_L, _environment$time), '\
            'Derivative(_L_type_Ca_channel_f_L_gate$f_L, _environment$time), '\
            'Derivative(_T_type_Ca_channel_d_T_gate$d_T, _environment$time), '\
            'Derivative(_T_type_Ca_channel_f_T_gate$f_T, _environment$time), '\
            'Derivative(_Ca_independent_transient_outward_K_current_r_gate$r, _environment$time), '\
            'Derivative(_Ca_independent_transient_outward_K_current_s1_gate$s1, _environment$time), '\
            'Derivative(_Ca_independent_transient_outward_K_current_s2_gate$s2, _environment$time), '\
            'Derivative(_Ca_independent_transient_outward_K_current_s3_gate$s3, _environment$time), '\
            'Derivative(_delayed_rectifier_K_current_z_gate$z, _environment$time), '\
            'Derivative(_delayed_rectifier_K_current_pa_gate$p_a, _environment$time), '\
            'Derivative(_delayed_rectifier_K_current_pi_gate$p_i, _environment$time), '\
            'Derivative(_intracellular_ion_concentrations$Na_i, _environment$time), '\
            'Derivative(_intracellular_ion_concentrations$Ca_i, _environment$time), '\
            'Derivative(_intracellular_ion_concentrations$K_i, _environment$time), '\
            'Derivative(_intracellular_Ca_buffering$O_C, _environment$time), '\
            'Derivative(_intracellular_Ca_buffering$O_TC, _environment$time), '\
            'Derivative(_intracellular_Ca_buffering$O_TMgC, _environment$time), '\
            'Derivative(_intracellular_Ca_buffering$O_TMgMg, _environment$time), '\
            'Derivative(_cleft_space_ion_concentrations$K_c, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$Ca_rel, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$Ca_up, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$O_Calse, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$F1, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$F2, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$F3, _environment$time)]'

    # also tested by model_units
    def test_get_equations_for(self, basic_model):
        """ Tests Model.get_equations_for() works correctly.
        Note: the basic model has 1 equation dsv11/dt = 1 which is not
        related to a symbol and so has no equations that can be retrieved
        by symbol
        """

        symbol_a = basic_model.get_symbol_by_cmeta_id("sv11")
        equation = basic_model.get_equations_for([symbol_a])
        assert len(equation) == 0

        symbol_t = basic_model.get_symbol_by_name("environment$time")

        equation1 = basic_model.get_equations_for([symbol_t])
        assert len(equation1) == 0

    # also tested by model_units
    def test_get_equations_for_1(self, aslanidi_model, OXMETA):
        """ Tests Model.get_equations_for() works correctly. """

        symbol_a = aslanidi_model.get_symbol_by_ontology_term(OXMETA, "membrane_capacitance")
        equation = aslanidi_model.get_equations_for([symbol_a])
        assert len(equation) == 1
        assert equation[0].lhs == symbol_a
        assert equation[0].rhs == 5e-5

    def test_get_value(self, aslanidi_model, OXMETA):
        """ Tests Model.get_value() works correctly. """

        symbol_a = aslanidi_model.get_symbol_by_ontology_term(OXMETA, "membrane_capacitance")
        assert aslanidi_model.get_value(symbol_a) == 5e-5

    #################################################################
    # tests for get_symbol_XXX functions

    def test_get_symbol_by_cmeta_id(self, basic_model):
        """ Tests Model.get_symbol_by_cmeta_id() works correctly. """

        sv11 = basic_model.get_symbol_by_cmeta_id('sv11')
        assert sv11.name == 'env_ode$sv1'
        assert sv11.units == 'mV'

    def test_get_symbol_by_cmeta_id_2(self, aslanidi_model):
        """ Tests Model.get_symbol_by_cmeta_id() works correctly. """

        variable = aslanidi_model.get_symbol_by_cmeta_id('testcmeta')
        assert variable.name == 'intracellular_ion_concentrations$Na_i'
        assert variable.units == 'millimolar'

    def test_get_symbol_by_name(self, basic_model):
        """ Tests Model.get_symbol_by_name() works correctly. """

        sv11 = basic_model.get_symbol_by_name('env_ode$sv1')
        assert sv11.units == 'mV'

    def test_get_symbol_by_ontology_term(self, aslanidi_model, OXMETA):
        """ Tests Model.get_symbol_by_ontology_term() works correctly. """

        symbol_a = aslanidi_model.get_symbol_by_ontology_term(OXMETA, 'membrane_capacitance')
        assert symbol_a.name == 'membrane$Cm'
        assert symbol_a.units == 'nanoF'

    def test_get_symbols_by_rdf(self, aslanidi_model, OXMETA):
        """ Tests Model.get_symbols_by_rdf() works correctly. """

        symbol_a = aslanidi_model.get_symbols_by_rdf(('http://biomodels.net/biology-qualifiers/', 'is'),
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
    # add_rdf

    #########################################################################
    # test add/set functions

    def test_add_equation(self, local_model):
        """ Tests the Model.add_equation method.
        """
        model = local_model
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

    def test_set_equation(self, local_model):
        """ Tests the Model.set_equation method.
        """
        model = local_model
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

    def test_set_equation2(self, local_hh_model, OXMETA):
        """ Tests replacing an equation in a model. """

        model = local_hh_model
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

    def test_add_number(self, local_model):
        """ Tests the Model.add_number method. """
        model = local_model
        number2 = model.add_number(2.0, 'mV')
        assert number2.is_Dummy

    def test_add_unit(self, local_model):
        """ Tests the Model.add_unit method. """
        model = local_model
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

    def test_add_variable(self, local_model):
        """ Tests the Model.add_variable() method. """
        model = local_model
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

    ###################################################################
    # this section is for other functions

    def test_find_symbols_and_derivatives(self, basic_model):
        """ Tests Model.find_symbols_and_derivatives function. """
        a = VariableDummy('a', 'second')
        b = VariableDummy('b', 'second')
        ex = sp.Add(a, b)
        syms = basic_model.find_symbols_and_derivatives([ex])
        assert len(syms) == 2

    def test_find_symbols_and_derivatives2(self, hh_model, OXMETA):
        """ Tests Model.find_symbols_and_derivatives() function. """

        # Test on single variable expressions
        t = hh_model.get_free_variable_symbol()
        syms = hh_model.find_symbols_and_derivatives([t])
        assert len(syms) == 1
        assert t in syms

        v = hh_model.get_symbol_by_ontology_term(OXMETA, "membrane_voltage")
        dvdt = sp.Derivative(v, t)
        syms = hh_model.find_symbols_and_derivatives([dvdt])
        assert len(syms) == 1
        assert sp.Derivative(v, t) in syms

        # Test on longer expressions
        x = sp.Float(1) + t * sp.sqrt(dvdt) - t
        syms = hh_model.find_symbols_and_derivatives([x])
        assert len(syms) == 2
        assert t in syms
        assert sp.Derivative(v, t) in syms

        # Test on multiple expressions
        y = sp.Float(2) + v
        syms = hh_model.find_symbols_and_derivatives([x, y])
        assert len(syms) == 3
        assert v in syms
        assert t in syms
        assert sp.Derivative(v, t) in syms

    def test_connect_variables(self, hh_model):
        """ Tests Model.connect_variables() function. """
        target = hh_model.get_symbol_by_name('sodium_channel$h')
        source = hh_model.get_symbol_by_name('sodium_channel_h_gate$h')
        assert target.assigned_to == source

        # check cannot assign already connected variable
        with pytest.raises(ValueError, match='Target already assigned'):
            hh_model.connect_variables('sodium_channel$h', 'sodium_channel_h_gate$h')

    def test_connect_variables2(self, local_hh_model):
        """ Tests Model.connect_variables() function. """
