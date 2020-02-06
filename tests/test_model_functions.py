import os

import pytest
import sympy as sp

from cellmlmanip import parser, units
from cellmlmanip.model import Model, VariableDummy
from cellmlmanip.model import FLOAT_PRECISION

from . import shared


class TestModelFunctions():
    """ Tests for all methods on Model class """

    ###############################################################
    # fixtures

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
        path = os.path.join(os.path.dirname(__file__), 'cellml_files', '4.algebraic_ode_model.cellml')
        with pytest.raises(ValueError, match='Equation LHS should be a derivative or variable'):
            parser.Parser(path).parse()

    #######################################################################
    # this section contains tests for each get_XXX function on Model

    def test_get_derived_quantities(self, basic_model, simple_ode_model):
        """ Tests Model.get_state_symbols() works correctly. """

        derived_quantities = basic_model.get_derived_quantities()
        assert len(derived_quantities) == 0

        derived_quantities = simple_ode_model.get_derived_quantities()
        assert str(derived_quantities) == \
            '[_single_ode_rhs_computed_var$a, _derived_from_state_var$dbl_sv1, _deriv_on_rhs$sv1_rate, '\
            '_circle_x_sibling$x2, _circle_y_implementation$rhs, _circle_sibling$local_complex_maths, '\
            '_time_units_conversion1$time, _deriv_on_rhs1a$sv1_rate, _time_units_conversion2$time, '\
            '_deriv_on_rhs2a$sv1_rate, _deriv_on_rhs1b$sv1_rate, _deriv_on_rhs2b$sv1_rate]'

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

    def test_get_equations_for(self):
        """
        Tests Model.get_equations_for().
        """

        m = Model('simplification')
        u = m.get_units('dimensionless')
        t = m.add_variable('t', u)
        y1 = m.add_variable('y1', u, initial_value=10)
        y2 = m.add_variable('y2', u, initial_value=20)
        y3 = m.add_variable('y3', u, initial_value=30)
        v1 = m.add_variable('v1', u)
        v2 = m.add_variable('v2', u)
        v3 = m.add_variable('v3', u)
        v4 = m.add_variable('v4', u)
        v5 = m.add_variable('v5', u)
        a1 = m.add_variable('a1', u)

        # dy1/dt = 1
        m.add_equation(sp.Eq(sp.Derivative(y1, t), m.add_number(1, u)))
        # dy2/dt = v1 --> Doesn't simplify, reference to v1 is maintained
        m.add_equation(sp.Eq(sp.Derivative(y2, t), v1))
        # dy3/dt = v2 * (2 + dy1/dt)
        m.add_equation(sp.Eq(sp.Derivative(y3, t), sp.Mul(v2, sp.Add(m.add_number(2, u), sp.Derivative(y1, t)))))
        # v1 = (5 - 5) * v3 --> Simplifies to 0
        m.add_equation(sp.Eq(v1, sp.Mul(sp.Add(m.add_number(5, u), m.add_number(-5, u)), v3)))
        # v2 = 23 + v4 --> Doesn't simplify, reference to v4 is maintained
        m.add_equation(sp.Eq(v2, sp.Add(m.add_number(23, u), v4)))
        # v3 = 2 / 3
        m.add_equation(sp.Eq(v3, sp.Mul(m.add_number(2, u), sp.Pow(m.add_number(3, u), sp.S.NegativeOne))))
        # v4 = -23
        m.add_equation(sp.Eq(v4, m.add_number(-23, u)))
        # v5 = v3 + v4
        m.add_equation(sp.Eq(v5, sp.Add(v3, v4)))
        # a1 = v5 + v2 + v1 + t
        m.add_equation(sp.Eq(a1, sp.Add(v5, v2, v1, t)))

        # Simplified equations
        e_v1 = sp.Eq(v1, sp.Float(0.0, FLOAT_PRECISION))
        e_v2 = sp.Eq(v2, sp.Add(v4, sp.Float(23., FLOAT_PRECISION)))
        e_v3 = sp.Eq(v3, sp.Float(2, FLOAT_PRECISION) / sp.Float(3, FLOAT_PRECISION))
        e_v4 = sp.Eq(v4, sp.Float(-23., FLOAT_PRECISION))
        e_v5 = sp.Eq(v5, sp.Add(v3, v4))
        e_a1 = sp.Eq(a1, sp.Add(v1, v2, v5, t))

        d_y1 = sp.Derivative(y1, t)
        d_y2 = sp.Derivative(y2, t)
        d_y3 = sp.Derivative(y3, t)

        e_y1 = sp.Eq(d_y1, sp.Float(1., FLOAT_PRECISION))
        e_y2 = sp.Eq(d_y2, v1)
        e_y3 = sp.Eq(d_y3, sp.Mul(v2, sp.Add(sp.Float(2., FLOAT_PRECISION), d_y1)))

        # v1 with simplification: [v1=0] (simplified)
        eqs = m.get_equations_for([v1])
        assert eqs[0] == e_v1
        assert len(eqs) == 1

        # v1 without simplification: [v3=2/3, v1=(5-5)*v3]
        eqs = m.get_equations_for([v1], strip_units=False)
        assert eqs[0] == m.graph.nodes[v3]['equation']
        assert eqs[1] == m.graph.nodes[v1]['equation']
        assert len(eqs) == 2

        # dy1/dt with simplification: [dy1/dt=1]
        eqs = m.get_equations_for([d_y1])
        assert eqs[0] == e_y1
        assert len(eqs) == 1

        # dy2/dt with simplification: [v1=0, dy2/dt=v1]
        eqs = m.get_equations_for([d_y2])
        assert eqs[0] == e_v1
        assert eqs[1] == e_y2
        assert len(eqs) == 2

        # dy2/dt without simplification: [v3=2/3, v1=(5-5)*v3, dy2/dt=v1]
        eqs = m.get_equations_for([d_y2], strip_units=False)
        assert eqs[0] == m.graph.nodes[v3]['equation']
        assert eqs[1] == m.graph.nodes[v1]['equation']
        assert eqs[2] == m.graph.nodes[d_y2]['equation']
        assert len(eqs) == 3

        # dy3/dt with simpification: [dy1/dt=1, v4=-23, v2=v4+23, dy2/dt=v2*(2+dy1/dt)]
        eqs = m.get_equations_for([d_y3])
        assert e_y3 in eqs
        assert e_y1 in eqs
        assert e_v2 in eqs
        assert e_v4 in eqs
        assert len(eqs) == 4

        # a1 with simplification: [v1=0, v3=2/3, v4=-23, v2=v4+23, v5=v3+v4, a1=v1+v2+v5]
        eqs = m.get_equations_for([a1])
        assert eqs[0] == e_v1
        assert eqs[1] == e_v3
        assert eqs[2] == e_v4
        assert eqs[3] == e_v2
        assert eqs[4] == e_v5
        assert eqs[5] == e_a1
        assert len(eqs) == 6

        # a1 with only one level of recursion
        eqs = m.get_equations_for([a1], recurse=False)
        assert eqs[0] == e_v1
        assert eqs[1] == e_v2
        assert eqs[2] == e_v5
        assert eqs[3] == e_a1
        assert len(eqs) == 4

        # Multiple input symbols: [d_y1=1, v1=0, d_y2=v1, v4=-23, v2=23+v4, d_y3=v2*(2+d_y1)]
        eqs = m.get_equations_for([d_y1, d_y2, d_y3])
        assert eqs[0] == e_y1
        assert eqs[1] == e_v1
        assert eqs[2] == e_y2
        assert eqs[3] == e_v4
        assert eqs[4] == e_v2
        assert eqs[5] == e_y3
        assert len(eqs) == 6

    def test_get_value(self, aslanidi_model):
        """ Tests Model.get_value() works correctly. """

        symbol_a = aslanidi_model.get_symbol_by_ontology_term(shared.OXMETA, "membrane_capacitance")
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

    def test_get_symbol_by_ontology_term(self, aslanidi_model):
        """ Tests Model.get_symbol_by_ontology_term() works correctly. """

        symbol_a = aslanidi_model.get_symbol_by_ontology_term(shared.OXMETA, 'membrane_capacitance')
        assert symbol_a.name == 'membrane$Cm'
        assert symbol_a.units == 'nanoF'

    def test_get_symbols_by_rdf(self, aslanidi_model):
        """ Tests Model.get_symbols_by_rdf() works correctly. """

        symbol_a = aslanidi_model.get_symbols_by_rdf(('http://biomodels.net/biology-qualifiers/', 'is'),
                                                     (shared.OXMETA, 'membrane_voltage'))
        assert len(symbol_a) == 1
        assert symbol_a[0].name == 'membrane$V'
        assert symbol_a[0].units == 'millivolt'

    ######################################################################
    # The functions listed for ontology/rdf are tested in test_rdf.py
    # Note these are functions that are not tested in this file as they are rdf related.
    #
    # In case of future changes this list was correct on 29 Nov 2019
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
        # but need to also add newvar1 = 2; newvar = 2 in order for the graph to resolve correctly
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

        # Now check some error cases
        with pytest.raises(ValueError, match='The variable newvar is defined twice'):
            # Redefining normal var
            model.add_equation(sp.Eq(symbol, 3.0))

        sv1 = model.get_symbol_by_cmeta_id('sv11')
        with pytest.raises(ValueError, match='The variable env_ode\\$sv1 is defined twice'):
            # Redefining state var
            model.add_equation(sp.Eq(sv1, 2.0))

        sv1_def = model._ode_definition_map[sv1]
        with pytest.raises(ValueError, match='The variable env_ode\\$sv1 is defined twice'):
            # Redefining with ODE
            model.add_equation(sp.Eq(sv1_def.lhs, 1.0))

        # But if we don't check for duplicates these are 'OK'
        model.add_equation(sp.Eq(symbol, 3.0), check_duplicates=False)
        model.add_equation(sp.Eq(sv1, 2.0), check_duplicates=False)
        model.add_equation(sp.Eq(sv1_def.lhs, 1.0), check_duplicates=False)

    def test_remove_equation(self, local_hh_model):
        """ Tests the Model.remove_equation method. """

        model = local_hh_model
        # Get model, assert that V is a state variable
        v = model.get_symbol_by_ontology_term(shared.OXMETA, 'membrane_voltage')
        # the issue here is that retrieving the variable uses the internal structure
        # which does not give the variable a type
        # to check type you need to use the graph
        # assert v.type == 'state'
        state_var = model.get_state_symbols()
        assert v in state_var

        # Now clamp it to -80mV
        t = model.get_symbol_by_ontology_term(shared.OXMETA, 'time')
        equation = model.graph.nodes[sp.Derivative(v, t)]['equation']
        model.remove_equation(equation)
        equation = sp.Eq(v, model.add_number(-80, v.units))
        model.add_equation(equation)

        # Check that V is no longer a state
        v = model.get_symbol_by_ontology_term(shared.OXMETA, 'membrane_voltage')
        state_var = model.get_state_symbols()
        assert v not in state_var

        # Now make V a state again
        dvdt_units = v.units / t.units
        model.remove_equation(equation)
        equation = sp.Eq(sp.Derivative(v, t), model.add_number(0, dvdt_units))
        model.add_equation(equation)

        # Check that V is a state again
        v = model.get_symbol_by_ontology_term(shared.OXMETA, 'membrane_voltage')
        state_var = model.get_state_symbols()
        assert v in state_var

        # Test removing non-existing equation
        equation = sp.Eq(sp.Derivative(v, t), model.add_number(5, dvdt_units))
        with pytest.raises(KeyError, match='Equation not found'):
            model.remove_equation(equation)

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
    # Unit related functionality

    def test_get_units(self):
        """ Tests Model.get_units(). """

        # Get predefined unit
        m = Model('test')
        m.get_units('volt')

        # Non-existent unit
        with pytest.raises(KeyError, match='Cannot find unit'):
            m.get_units('towel')

    def test_units(self, simple_units_model):
        """ Tests units read and calculated from a model. """
        symbol_a = simple_units_model.get_symbol_by_cmeta_id("a")
        equation = simple_units_model.get_equations_for([symbol_a], strip_units=False)
        assert simple_units_model.units.summarise_units(equation[0].lhs) == 'ms'
        assert simple_units_model.units.summarise_units(equation[0].rhs) == 'ms'

        symbol_b = simple_units_model.get_symbol_by_cmeta_id("b")
        equation = simple_units_model.get_equations_for([symbol_b])
        assert simple_units_model.units.summarise_units(equation[1].lhs) == 'per_ms'
        assert simple_units_model.units.summarise_units(equation[1].rhs) == '1 / ms'
        assert simple_units_model.units.is_unit_equal(simple_units_model.units.summarise_units(equation[1].lhs),
                                                      simple_units_model.units.summarise_units(equation[1].rhs))

    def test_bad_units(self, bad_units_model):
        """ Tests units read and calculated from an inconsistent model. """
        symbol_a = bad_units_model.get_symbol_by_cmeta_id("a")
        symbol_b = bad_units_model.get_symbol_by_cmeta_id("b")
        equation = bad_units_model.get_equations_for([symbol_b], strip_units=False)
        assert len(equation) == 2
        assert equation[0].lhs == symbol_a
        assert bad_units_model.units.summarise_units(equation[0].lhs) == 'ms'
        with pytest.raises(units.UnitError):
            # cellml file states a (ms) = 1 (ms) + 1 (second)
            bad_units_model.units.summarise_units(equation[0].rhs)

        assert equation[1].lhs == symbol_b
        with pytest.raises(units.UnitError):
            # cellml file states b (per_ms) = power(a (ms), 1 (second))
            bad_units_model.units.summarise_units(equation[1].rhs)

    ###################################################################
    # this section is for other functions

    def test_find_symbols_and_derivatives(self, basic_model):
        """ Tests Model.find_symbols_and_derivatives function. """
        a = VariableDummy('a', 'second')
        b = VariableDummy('b', 'second')
        ex = sp.Add(a, b)
        syms = basic_model.find_symbols_and_derivatives([ex])
        assert len(syms) == 2

    def test_find_symbols_and_derivatives2(self, hh_model):
        """ Tests Model.find_symbols_and_derivatives() function. """

        # Test on single variable expressions
        t = hh_model.get_free_variable_symbol()
        syms = hh_model.find_symbols_and_derivatives([t])
        assert len(syms) == 1
        assert t in syms

        v = hh_model.get_symbol_by_ontology_term(shared.OXMETA, "membrane_voltage")
        dvdt = sp.Derivative(v, t)
        syms = hh_model.find_symbols_and_derivatives([dvdt])
        assert len(syms) == 1
        assert sp.Derivative(v, t) in syms

        # Test on longer expressions
        x = sp.Float(1, FLOAT_PRECISION) + t * sp.sqrt(dvdt) - t
        syms = hh_model.find_symbols_and_derivatives([x])
        assert len(syms) == 2
        assert t in syms
        assert sp.Derivative(v, t) in syms

        # Test on multiple expressions
        y = sp.Float(2, FLOAT_PRECISION) + v
        syms = hh_model.find_symbols_and_derivatives([x, y])
        assert len(syms) == 3
        assert v in syms
        assert t in syms
        assert sp.Derivative(v, t) in syms

    def test_connect_variables(self, local_hh_model):
        """ Tests Model.connect_variables() function. """
        target = local_hh_model.get_symbol_by_name('sodium_channel$h')
        source = local_hh_model.get_symbol_by_name('sodium_channel_h_gate$h')
        assert target.assigned_to == source

        # check cannot assign already connected variable
        with pytest.raises(ValueError, match='Target already assigned'):
            local_hh_model.connect_variables('sodium_channel$h', 'sodium_channel_h_gate$h')

        # add and connect a new variable with same units
        num_eqns = len(local_hh_model.equations)
        local_hh_model.add_variable(name='newvar', units='dimensionless', public_interface='in')
        new_target = local_hh_model.get_symbol_by_name('newvar')
        local_hh_model.connect_variables('sodium_channel_h_gate$h', 'newvar')
        assert new_target.assigned_to == source
        assert len(local_hh_model.equations) == num_eqns

        # Can't connect a variable already defined by an equation to another source
        newvar2 = local_hh_model.add_variable(name='newvar2', units='dimensionless', public_interface='in')
        local_hh_model.add_equation(sp.Eq(newvar2, 1.0))
        with pytest.raises(ValueError, match='Multiple definitions for newvar2'):
            local_hh_model.connect_variables('sodium_channel$E_Na', 'newvar2')

    def test_connect_variable2(self, local_hh_model):
        """ Tests Model.connect_variables() function. """
        num_eqns = len(local_hh_model.equations)
        # add and connect a variable that requires unit conversion
        local_hh_model.add_variable(name='newvar', units='volt', public_interface='in')
        new_target = local_hh_model.get_symbol_by_name('newvar')
        source = local_hh_model.get_symbol_by_name('leakage_current$E_L')
        local_hh_model.connect_variables('leakage_current$E_L', 'newvar')
        assert new_target.assigned_to == new_target
        assert len(local_hh_model.equations) == num_eqns + 1
        new_eqn = local_hh_model.equations[-1]
        assert new_eqn.is_Equality
        assert new_eqn.lhs == new_target
        assert new_eqn.rhs.is_Mul
        assert new_eqn.rhs.args[0].value == 0.001
        assert new_eqn.rhs.args[1] == source
