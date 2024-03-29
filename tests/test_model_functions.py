import os

import pytest
import sympy as sp

from cellmlmanip import parser, units
from cellmlmanip.model import FLOAT_PRECISION, Model, Variable

from . import shared


class TestModelFunctions():
    """
    Tests for most methods in Model.

    RDF-related methods are tested in test_rdf.py
    Unit-related methods in test_unit_conversion.py
    """

    ###############################################################
    # fixtures

    @pytest.fixture
    def local_model(scope='function'):
        """ Fixture to load a local copy of the basic_ode model that may get modified. """
        return shared.load_model('basic_ode')

    @pytest.fixture
    def local_hh_model(scope='function'):
        """ Fixture to load a local copy of the hodgkin_huxley_squid_axon_model_1952_modified
        model that may get modified. """
        return shared.load_model('hodgkin_huxley_squid_axon_model_1952_modified')

    ##########################################################
    # check equation graph property

    def test_graph_property(self, basic_model):
        """ Tests that the graph property for Model has been constructed correctly. """

        graph1 = basic_model.graph
        assert len(graph1.nodes) == 3
        names = ['env_ode$sv1', 'environment$time']
        for v in graph1:
            if not v.is_Derivative:
                assert v.name in names
            else:
                for a in v._args:
                    if a.is_Dummy:
                        assert a.name in names
                    else:
                        for b in a._args:
                            if b.is_Dummy:
                                assert b.name in names

    def test_graph_for_dae(self):
        """ Checks if writing a DAE in a model raises an exceptions. """
        path = os.path.join(os.path.dirname(__file__), 'cellml_files', '4.algebraic_ode_model.cellml')
        with pytest.raises(ValueError, match='Equation LHS should be a derivative or variable'):
            parser.Parser(path).parse()

    #######################################################################
    # this section contains tests for each get_XXX function on Model

    def test_get_derived_quantities(self, basic_model, simple_ode_model):
        """ Tests Model.get_derived_quantities(). """

        derived_quantities = basic_model.get_derived_quantities()
        unsorted_derived_quantities = basic_model.get_derived_quantities(sort=False)
        assert len(derived_quantities) == len(unsorted_derived_quantities) == 0

        derived_quantities = simple_ode_model.get_derived_quantities()
        unsorted_derived_quantities = simple_ode_model.get_derived_quantities(sort=False)
        assert set(derived_quantities) == set(unsorted_derived_quantities)
        assert str(derived_quantities) == (
            '['
            # '_single_ode_rhs_computed_var$a, '
            '_derived_from_state_var$dbl_sv1, '
            '_deriv_on_rhs$sv1_rate, '
            '_circle_x_sibling$x2, '
            '_circle_y_implementation$rhs, '
            '_circle_sibling$local_complex_maths, '
            '_time_units_conversion1$time, '
            '_deriv_on_rhs1a$sv1_rate, '
            '_time_units_conversion2$time, '
            '_deriv_on_rhs2a$sv1_rate, '
            '_deriv_on_rhs1b$sv1_rate, '
            '_deriv_on_rhs2b$sv1_rate]'
        )

    def test_get_display_name(self, simple_ode_model):
        """ Tests Model.get_display_name(var). """
        OXMETA = 'https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#'

        # For getting a display name there are 4 possibilities:
        # Has an oxmeta-annotation
        var = simple_ode_model.get_variable_by_name('single_independent_ode$sv1')
        assert var.cmeta_id == 'sv11'
        assert simple_ode_model.get_display_name(var, ontology=OXMETA) == 'sodium_reversal_potential'

        # Has an annotation in another ontology
        var = simple_ode_model.get_variable_by_name('single_ode_rhs_const_var$a')
        assert var.cmeta_id == 'a1'
        assert simple_ode_model.get_display_name(var, ontology=OXMETA) == var.cmeta_id
        assert simple_ode_model.get_display_name(var) == 'parameter_a1'
        assert simple_ode_model.get_display_name(var, ontology='urn:test-ns#') == 'parameter_a1'

        # Has no annotation but does have a cmeta:id
        var = simple_ode_model.get_variable_by_name('single_ode_rhs_computed_var$a')
        assert var.cmeta_id == 'a2'
        assert simple_ode_model.get_display_name(var, OXMETA) == var.cmeta_id
        assert simple_ode_model.get_display_name(var) == var.cmeta_id

        # Has no cmeta:id
        var = simple_ode_model.get_variable_by_name('single_ode_rhs_const_var$time')
        assert var.cmeta_id is None
        assert simple_ode_model.get_display_name(var) == var.name.replace('$', '__')

        # Check excluded tags work
        var = simple_ode_model.get_variable_by_name('single_independent_ode$sv1')
        assert var.cmeta_id == 'sv11'
        assert simple_ode_model.get_display_name(var, OXMETA, set({'sodium_reversal_potential'})) == 'sv11'

    def test_get_state_variables(self, basic_model):
        """ Tests Model.get_state_variables() works on a simple model. """

        states = basic_model.get_state_variables()
        unsorted_states = basic_model.get_state_variables(sort=False)
        assert states == unsorted_states
        assert len(states) == 1
        assert states[0].name == 'env_ode$sv1'

    def test_get_state_variables_2(self, aslanidi_model):
        """ Tests Model.get_state_variables() works on a complex model. """
        states = aslanidi_model.get_state_variables()
        unsorted_states = aslanidi_model.get_state_variables(sort=False)
        assert set(states) == set(unsorted_states)
        assert len(states) == 29
        assert str(states) == \
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
    def test_get_free_variable(self, basic_model, aslanidi_model):
        """ Tests Model.get_free_variable() works correctly. """

        var = basic_model.get_free_variable()
        assert var.name == 'environment$time'

        var = aslanidi_model.get_free_variable()
        assert var.name == 'environment$time'

    def test_get_derivatives(self, basic_model):
        """ Tests Model.get_derivatives() works correctly. """

        derivs = basic_model.get_derivatives()
        unsorted_derivs = basic_model.get_derivatives(sort=False)
        assert set(derivs) == set(unsorted_derivs)
        assert len(derivs) == 1
        deriv = derivs[0]
        assert deriv.is_Derivative
        assert len(deriv.variables) == 1
        assert deriv.variables[0].is_Dummy
        assert deriv.variables[0].name == 'environment$time'

    def test_get_derivatives_2(self, aslanidi_model):
        """ Tests Model.get_derivatives() works correctly on a more complicated model. """

        derivs = aslanidi_model.get_derivatives()
        unsorted_derivs = aslanidi_model.get_derivatives(sort=False)
        assert set(derivs) == set(unsorted_derivs)
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

    def test_unsupported_unit(self, aslanidi_model):
        with pytest.raises(KeyError):
            aslanidi_model.units.get_unit('celsius')

    def test_get_equations_for(self):
        """
        Tests Model.get_equations_for().
        """

        m = Model('simplification')
        u = m.units.get_unit('dimensionless')
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
        m.add_equation(sp.Eq(sp.Derivative(y1, t), m.create_quantity(1, u)))
        # dy2/dt = v1 --> Doesn't simplify, reference to v1 is maintained
        m.add_equation(sp.Eq(sp.Derivative(y2, t), v1))
        # dy3/dt = v2 * (2 + dy1/dt)
        m.add_equation(sp.Eq(sp.Derivative(y3, t), sp.Mul(v2, sp.Add(m.create_quantity(2, u), sp.Derivative(y1, t)))))
        # v1 = (5 - 5) * v3 --> Simplifies to 0
        m.add_equation(sp.Eq(v1, sp.Mul(sp.Add(m.create_quantity(5, u), m.create_quantity(-5, u)), v3)))
        # v2 = 23 + v4 --> Doesn't simplify, reference to v4 is maintained
        m.add_equation(sp.Eq(v2, sp.Add(m.create_quantity(23, u), v4)))
        # v3 = 2 / 3
        m.add_equation(sp.Eq(v3, sp.Mul(m.create_quantity(2, u), sp.Pow(m.create_quantity(3, u), sp.S.NegativeOne))))
        # v4 = -23
        m.add_equation(sp.Eq(v4, m.create_quantity(-23, u)))
        # v5 = v3 + v4
        m.add_equation(sp.Eq(v5, sp.Add(v3, v4)))
        # a1 = v5 + v2 + v1 + t
        m.add_equation(sp.Eq(a1, sp.Add(v5, v2, v1, t)))

        # Simplified equations
        e_v1 = sp.Eq(v1, sp.Number(0))
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

    def test_get_definition(self, simple_ode_model):
        # Simple equation
        a = simple_ode_model.get_variable_by_cmeta_id('a2')
        defn = simple_ode_model.get_definition(a)
        assert str(defn) == 'Eq(_single_ode_rhs_computed_var$a, _-1)'

        # ODE definition
        state_var = simple_ode_model.get_variable_by_cmeta_id('sv11')
        defn = simple_ode_model.get_definition(state_var)
        assert str(defn) == 'Eq(Derivative(_single_independent_ode$sv1, _environment$time), _1)'

        # No defining equation
        time = simple_ode_model.get_variable_by_cmeta_id('time')
        defn = simple_ode_model.get_definition(time)
        assert defn is None

    def test_get_value(self, local_hh_model):
        """ Tests Model.get_value() works correctly. """
        model = local_hh_model

        # Test getting a constant
        C = model.get_variable_by_ontology_term((shared.OXMETA, 'membrane_capacitance'))
        assert model.get_value(C) == 1.0

        # Test getting a state
        V = model.get_variable_by_ontology_term((shared.OXMETA, 'membrane_voltage'))
        assert model.get_value(V) == V.initial_value

        # Test getting an intermediary value
        INa = model.get_variable_by_ontology_term((shared.OXMETA, 'membrane_fast_sodium_current'))
        assert model.get_value(INa) == pytest.approx(-1.03500000000000014e+00)

        # Test getting the time variable
        t = model.get_free_variable()
        assert model.get_value(t) == 0

        # Test getting a value that depends on time
        model.remove_equation(model.get_definition(C))
        model.add_equation(sp.Eq(C, 4 * (0.5 + t)))
        assert model.get_value(C) == 2

        # Test getting a value in a model without states
        for state in model.get_state_variables():
            model.remove_equation(model.get_definition(state))
            model.add_equation(sp.Eq(state, model.create_quantity(state.initial_value, state.units)))
        assert model.get_value(INa) == pytest.approx(-1.03500000000000014e+00)

        # Test getting an undefined variable
        model.remove_equation(model.get_definition(C))
        with pytest.raises(ValueError, match='No definition'):
            model.get_value(C)

    #################################################################
    # tests for get_variable_XXX functions

    def test_get_variable_by_cmeta_id(self, basic_model):
        """ Tests Model.get_variable_by_cmeta_id() works correctly. """
        sv11 = basic_model.get_variable_by_cmeta_id('sv11')
        assert str(sv11.rdf_identity) == '#sv11'
        assert sv11.name == 'env_ode$sv1'
        assert sv11.units == basic_model.units.get_unit('mV')

        term = sv11.rdf_identity
        assert basic_model.get_variable_by_cmeta_id(term) is sv11

    def test_get_variable_by_name(self, basic_model):
        """ Tests Model.get_variable_by_name() works correctly. """

        sv11 = basic_model.get_variable_by_name('env_ode$sv1')
        assert sv11.units == basic_model.units.get_unit('mV')

    def test_get_variable_by_ontology_term(self, aslanidi_model):
        """ Tests Model.get_variable_by_ontology_term() works correctly. """

        symbol_a = aslanidi_model.get_variable_by_ontology_term((shared.OXMETA, 'membrane_capacitance'))
        assert symbol_a.name == 'membrane$Cm'
        assert symbol_a.units == aslanidi_model.units.get_unit('nanoF')

    def test_get_variables_by_rdf(self, aslanidi_model):
        """ Tests Model.get_variables_by_rdf() works correctly. """

        predicate = ('http://biomodels.net/biology-qualifiers/', 'is')
        object_ = (shared.OXMETA, 'membrane_voltage')

        symbol_a = aslanidi_model.get_variables_by_rdf(predicate, object_)
        unsorted_symbol_a = aslanidi_model.get_variables_by_rdf(predicate, object_, sort=False)
        assert set(symbol_a) == set(unsorted_symbol_a)
        assert len(symbol_a) == 1
        assert symbol_a[0].name == 'membrane$V'
        assert symbol_a[0].units == aslanidi_model.units.get_unit('millivolt')

    ######################################################################
    # The functions listed for ontology/rdf are tested in test_rdf.py
    # Note these are functions that are not tested in this file as they are rdf related.
    #
    # In case of future changes this list was correct on 29 Nov 2019
    #
    # get_ontology_terms_by_variable()
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
        symbol = model.get_variable_by_name('newvar')
        model.add_variable(name='newvar1', units='mV')
        symbol1 = model.get_variable_by_name('newvar1')
        model.add_variable(name='newvar2', units='mV')
        symbol2 = model.get_variable_by_name('newvar2')
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

        sv1 = model.get_variable_by_cmeta_id('sv11')
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

        # Only first-order derivatives allowed
        with pytest.raises(ValueError, match='Only first order derivatives'):
            model.add_equation(sp.Eq(sp.Derivative(symbol, symbol1, symbol2), 1.0))
        with pytest.raises(ValueError, match='Only first order derivatives'):
            model.add_equation(sp.Eq(sp.Derivative(symbol, symbol1, symbol1), 1.0))

    def test_remove_equation(self, local_hh_model):
        """ Tests the Model.remove_equation method. """

        model = local_hh_model
        # Get model, assert that V is a state variable
        v = model.get_variable_by_ontology_term((shared.OXMETA, 'membrane_voltage'))
        # the issue here is that retrieving the variable uses the internal structure
        # which does not give the variable a type
        # to check type you need to use the graph
        # assert v.type == 'state'
        state_var = model.get_state_variables()
        assert v in state_var

        # Now clamp it to -80mV
        t = model.get_variable_by_ontology_term((shared.OXMETA, 'time'))
        equation = model.graph.nodes[sp.Derivative(v, t)]['equation']
        model.remove_equation(equation)
        equation = sp.Eq(v, model.create_quantity(-80, v.units))
        model.add_equation(equation)

        # Check that V is no longer a state
        v = model.get_variable_by_ontology_term((shared.OXMETA, 'membrane_voltage'))
        state_var = model.get_state_variables()
        assert v not in state_var

        # Now make V a state again
        dvdt_units = v.units / t.units
        model.remove_equation(equation)
        equation = sp.Eq(sp.Derivative(v, t), model.create_quantity(0, dvdt_units))
        model.add_equation(equation)

        # Check that V is a state again
        v = model.get_variable_by_ontology_term((shared.OXMETA, 'membrane_voltage'))
        state_var = model.get_state_variables()
        assert v in state_var

        # Test removing non-existing equation
        equation = sp.Eq(sp.Derivative(v, t), model.create_quantity(5, dvdt_units))
        with pytest.raises(KeyError, match='Equation not found'):
            model.remove_equation(equation)

    def test_create_quantity(self, local_model):
        """ Tests the Model.create_quantity method. """
        model = local_model
        number2 = model.create_quantity(2.0, 'mV')
        assert number2.is_Dummy

    def test_add_variable(self, local_model):
        """ Tests the Model.add_variable() method. """
        model = local_model
        assert len(model.variables()) == 3
        with pytest.raises(KeyError):
            model.get_variable_by_name('newvar') is None

        newvar = model.add_variable('newvar', 'mV')
        assert newvar.model is model
        assert newvar.is_real
        assert len(model.variables()) == 4
        assert model.get_variable_by_name('newvar')

        # Variable can't be added twice
        model.add_variable('varvar1', 'mV')
        model.add_variable('varvar2', 'mV')
        assert len(model.variables()) == 6
        with pytest.raises(ValueError, match='already exists'):
            model.add_variable('varvar1', 'mV')

        # Same cmeta id can't be used twice
        model.add_variable('v1', 'mV', cmeta_id='v1')
        with pytest.raises(ValueError, match='cmeta id'):
            model.add_variable('v2', 'mV', cmeta_id='v1')

    def test_remove_variable(self, local_hh_model):
        """Test the Model.remove_variable() method."""
        model = local_hh_model
        # Removing a 'normal' variable
        i_stim = model.get_variable_by_cmeta_id('membrane_stimulus_current')
        i_stim_defn = model.get_definition(i_stim)
        assert i_stim_defn is not None
        i_stim_rdf = list(model.get_rdf_annotations(i_stim.rdf_identity))
        assert len(i_stim_rdf) == 1

        model.remove_variable(i_stim)

        assert i_stim.model is None
        with pytest.raises(KeyError):
            model.get_variable_by_name(i_stim.name)
        assert model.get_definition(i_stim) is None
        assert i_stim_defn not in model.equations
        assert len(list(model.get_rdf_annotations(i_stim.rdf_identity))) == 0
        assert i_stim_rdf[0] not in model.rdf

        # Removing a state variable
        v = model.get_variable_by_cmeta_id('membrane_voltage')
        v_defn = model.get_definition(v)
        assert v_defn is not None
        v_rdf = list(model.get_rdf_annotations(v.rdf_identity))
        assert len(v_rdf) == 1

        model.remove_variable(v)

        assert v.model is None
        with pytest.raises(KeyError):
            model.get_variable_by_name(v.name)
        assert model.get_definition(v) is None
        assert v_defn not in model.equations
        assert len(list(model.get_rdf_annotations(v.rdf_identity))) == 0
        assert v_rdf[0] not in model.rdf

        # Remove a variable with no definition
        t = model.get_variable_by_cmeta_id('time')
        assert model.get_definition(t) is None

        model.remove_variable(t)

        assert t.model is None
        with pytest.raises(KeyError):
            model.get_variable_by_name(t.name)

        # Remove a variable with no cmeta_id
        var = model.get_variable_by_name('membrane$E_R')
        model.remove_variable(var)

        assert var.model is None
        with pytest.raises(KeyError):
            model.get_variable_by_name(var.name)

    ###################################################################
    # Unit related functionality

    def test_units(self, simple_units_model):
        """ Tests units read and calculated from a model. """
        symbol_a = simple_units_model.get_variable_by_cmeta_id("a")
        equation = simple_units_model.get_equations_for([symbol_a], strip_units=False)
        assert simple_units_model.units.evaluate_units(equation[0].lhs) == simple_units_model.units.get_unit('ms')
        assert simple_units_model.units.evaluate_units(equation[0].rhs) == simple_units_model.units.get_unit('ms')

        symbol_b = simple_units_model.get_variable_by_cmeta_id("b")
        equation = simple_units_model.get_equations_for([symbol_b])
        assert simple_units_model.units.evaluate_units(equation[1].lhs) == simple_units_model.units.get_unit('per_ms')
        assert simple_units_model.units.evaluate_units(equation[1].rhs) == 1 / simple_units_model.units.get_unit('ms')
        assert simple_units_model.units.is_equivalent(
            simple_units_model.units.evaluate_units(equation[1].lhs),
            simple_units_model.units.evaluate_units(equation[1].rhs))

    def test_bad_units(self, bad_units_model):
        """ Tests units read and calculated from an inconsistent model. """
        symbol_a = bad_units_model.get_variable_by_cmeta_id("a")
        symbol_b = bad_units_model.get_variable_by_cmeta_id("b")
        equation = bad_units_model.get_equations_for([symbol_b], strip_units=False)
        assert len(equation) == 2
        assert equation[0].lhs == symbol_a
        assert bad_units_model.units.evaluate_units(equation[0].lhs) == bad_units_model.units.get_unit('ms')
        with pytest.raises(units.UnitError):
            # cellml file states a (ms) = 1 (ms) + 1 (second)
            bad_units_model.units.evaluate_units(equation[0].rhs)

        assert equation[1].lhs == symbol_b
        with pytest.raises(units.UnitError):
            # cellml file states b (per_ms) = power(a (ms), 1 (second))
            bad_units_model.units.evaluate_units(equation[1].rhs)

    ###################################################################
    # this section is for other functions

    def test_find_variables_and_derivatives(self, basic_model):
        """ Tests Model.find_variables_and_derivatives() on a simple model. """
        a = Variable('a', 'second')
        b = Variable('b', 'second')
        ex = sp.Add(a, b)
        syms = basic_model.find_variables_and_derivatives([ex])
        assert len(syms) == 2

    def test_find_variables_and_derivatives_2(self, hh_model):
        """ Tests Model.find_variables_and_derivatives() on a more complicated model. """

        # Test on single variable expressions
        t = hh_model.get_free_variable()
        syms = hh_model.find_variables_and_derivatives([t])
        assert len(syms) == 1
        assert t in syms

        v = hh_model.get_variable_by_ontology_term((shared.OXMETA, "membrane_voltage"))
        dvdt = sp.Derivative(v, t)
        syms = hh_model.find_variables_and_derivatives([dvdt])
        assert len(syms) == 1
        assert sp.Derivative(v, t) in syms

        # Test on longer expressions
        x = sp.Float(1, FLOAT_PRECISION) + t * sp.sqrt(dvdt) - t
        syms = hh_model.find_variables_and_derivatives([x])
        assert len(syms) == 2
        assert t in syms
        assert sp.Derivative(v, t) in syms

        # Test on multiple expressions
        y = sp.Float(2, FLOAT_PRECISION) + v
        syms = hh_model.find_variables_and_derivatives([x, y])
        assert len(syms) == 3
        assert v in syms
        assert t in syms
        assert sp.Derivative(v, t) in syms

    def test_variable_classification(self, aslanidi_model):
        """ Tests Model.is_state() and Model.is_constant(). """

        # State variable
        v = aslanidi_model.get_variable_by_ontology_term((shared.OXMETA, 'membrane_voltage'))
        assert aslanidi_model.is_state(v)
        assert not aslanidi_model.is_constant(v)

        # Intermediary variable
        i = aslanidi_model.get_variable_by_ontology_term((shared.OXMETA, 'membrane_fast_sodium_current'))
        assert not aslanidi_model.is_state(i)
        assert not aslanidi_model.is_constant(i)

        # Free variable
        t = aslanidi_model.get_variable_by_ontology_term((shared.OXMETA, 'time'))
        assert not aslanidi_model.is_state(t)
        assert not aslanidi_model.is_constant(t)

        # Constant
        c = aslanidi_model.get_variable_by_ontology_term((shared.OXMETA, 'membrane_capacitance'))
        assert not aslanidi_model.is_state(c)
        assert aslanidi_model.is_constant(c)
