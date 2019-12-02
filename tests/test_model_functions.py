import os

import pytest
import sympy as sp

from cellmlmanip import parser
from cellmlmanip.model import Model, VariableDummy

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

    @pytest.fixture
    def local_model_with_simplification(scope='function'):
        """ Fixture to load a model where simplification can be applied by sympy. """

        m = Model('simplification')
        u = 'dimensionless'
        t = m.add_variable('t', u)
        y1 = m.add_variable('y1', u, initial_value=10)
        y2 = m.add_variable('y2', u, initial_value=20)
        y3 = m.add_variable('y3', u, initial_value=30)
        v1 = m.add_variable('v1', u)
        v2 = m.add_variable('v2', u)
        v3 = m.add_variable('v3', u)
        v4 = m.add_variable('v4', u)
        v5 = m.add_variable('v5', u)

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

        return m

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
    def test_get_initial_value(self, aslanidi_model):
        """ Tests Model.get_initial_value() works correctly. """

        membrane_voltage = aslanidi_model.get_symbol_by_ontology_term(shared.OXMETA, "membrane_voltage")
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

    def test_get_equations_for(self, hh_model):
        """
        Tests Model.get_equations_for(), using ``strip_units=True`` (the default). """

        # Test get_equations_for with topgraphical lexicographical ordering

        # Get ordered equations
        membrane_fast_sodium_current = hh_model.get_symbol_by_ontology_term(shared.OXMETA,
                                                                            'membrane_fast_sodium_current')
        equations = hh_model.get_equations_for([membrane_fast_sodium_current])
        top_level_equations = hh_model.get_equations_for([membrane_fast_sodium_current], recurse=False)

        # There should be 4 in this model
        assert len(equations) == 4

        # There should be 3 top level (non recursed) in this model
        assert len(top_level_equations) == 3

        # Expected equations
        ER = sp.Eq(sp.Dummy('membrane$E_R'), sp.numbers.Float(-75.0))
        ENa = sp.Eq(sp.Dummy('sodium_channel$E_Na'),
                    sp.add.Add(sp.Dummy('membrane$E_R'), sp.numbers.Float(115.0)))
        gNa = sp.Eq(sp.Dummy('sodium_channel$g_Na'), sp.numbers.Float(120.0))
        iNa = sp.Eq(sp.Dummy('sodium_channel$i_Na'),
                    sp.Dummy('sodium_channel_m_gate$m') ** 3.0 * sp.Dummy('sodium_channel$g_Na') *
                    sp.Dummy('sodium_channel_h_gate$h') * (sp.Dummy('membrane$V') -
                    sp.Dummy('sodium_channel$E_Na')))

        # Get order as strings, for easy comparison
        ER, ENa, gNa, iNa = str(ER), str(ENa), str(gNa), str(iNa)
        expected_order = [ER, ENa, gNa, iNa]

        # Check equations against expected equations
        equations = [str(eq) for eq in equations]
        assert equations == expected_order

        # Check topologically (but not lexicographically) ordered equations
        unordered_equations = hh_model.get_equations_for([membrane_fast_sodium_current], False)
        unordered_equations = [str(eq) for eq in unordered_equations]

        # Each equation should be both in the ordered and unordered equations
        assert set(unordered_equations) == set(equations)

        # ER should come before ENa
        assert unordered_equations.index(ER) < unordered_equations.index(ENa)

        # ENa and gNa should come before iNa
        assert unordered_equations.index(ENa) < unordered_equations.index(iNa)
        assert unordered_equations.index(gNa) < unordered_equations.index(iNa)

    def test_get_equations_for_with_simplification(self, local_model_with_simplification):
        """
        Tests Model.get_equations_for() in a situation where sympy can make simplifications.
        """
        m = local_model_with_simplification

        # Get equations for v1: [v1=0] (simplified)
        v1 = m.get_symbol_by_name('v1')
        e_v1 = sp.Eq(v1, sp.Number(0))
        eqs = m.get_equations_for([v1])
        assert e_v1 in eqs
        assert len(eqs) == 1

        # Get equations for v2: [v2=v4+23, v4=-23]
        v4 = m.get_symbol_by_name('v4')
        e_v4 = sp.Eq(v4, sp.Number(-23))
        v2 = m.get_symbol_by_name('v2')
        e_v2 = sp.Eq(v2, sp.Add(v4, sp.Number(23)))
        eqs = m.get_equations_for([v2])
        assert e_v4 in eqs
        assert e_v2 in eqs
        assert len(eqs) == 2

        # Get equations for v3: [v3=0.67]
        v3 = m.get_symbol_by_name('v3')
        e_v3 = sp.Eq(v3, sp.Number(2 / 3))
        eqs = m.get_equations_for([v3])
        assert e_v3 in eqs
        assert len(eqs) == 1

        # Get equations for v4: [v4=-23]
        eqs = m.get_equations_for([v4])
        assert e_v4 in eqs
        assert len(eqs) == 1

        # Get equations for v5: [v5=v3+v4, v3=0.67, v4=-23]
        v5 = m.get_symbol_by_name('v5')
        e_v5 = sp.Eq(v5, sp.Add(v3, v4))
        eqs = m.get_equations_for([v5])
        assert e_v3 in eqs
        assert e_v4 in eqs
        assert e_v5 in eqs
        assert len(eqs) == 3

        # Get equations for dy1/dt: [dy1/dt=1]
        t = m.get_symbol_by_name('t')
        y1 = m.get_symbol_by_name('y1')
        d_y1 = sp.Derivative(y1, t)
        e_y1 = sp.Eq(d_y1, sp.Number(1))
        eqs = m.get_equations_for([d_y1])
        assert e_y1 in eqs
        assert len(eqs) == 1

        # Get equations for dy2/dt: [dy2/dt=v1, v1=0]
        y2 = m.get_symbol_by_name('y2')
        d_y2 = sp.Derivative(y2, t)
        e_y2 = sp.Eq(d_y2, v1)
        eqs = m.get_equations_for([d_y2])
        assert e_v1 in eqs
        assert e_y2 in eqs
        assert len(eqs) == 2

        # Get equations for dy3/dt: [dy2/dt=v2*(2+dy1/dt), dy1/dt=1, v2=v4+23, v4=-23]
        y3 = m.get_symbol_by_name('y3')
        d_y3 = sp.Derivative(y3, t)
        e_y3 = sp.Eq(d_y3, sp.Mul(v2, sp.Add(sp.Number(2), d_y1)))
        eqs = m.get_equations_for([d_y3])
        assert e_y3 in eqs
        assert e_y1 in eqs
        assert e_v2 in eqs
        assert e_v4 in eqs
        assert len(eqs) == 4

    def test_get_equations_for_with_dummies(self, hh_model):
        """
        Tests Model.get_equations_for() using ``strip_units=False``.
        """

        # Get ordered equations
        ina = hh_model.get_symbol_by_ontology_term(shared.OXMETA, 'membrane_fast_sodium_current')
        equations = hh_model.get_equations_for([ina], recurse=False, strip_units=False)

        # Test dummies are preserved
        for eq in equations:
            if eq.lhs.name == 'sodium_channel$g_Na':
                assert isinstance(eq.rhs, sp.Dummy)
                break

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

    def test_set_equation2(self, local_hh_model):
        """ Tests replacing an equation in a model. """

        model = local_hh_model
        # Get model, assert that V is a state variable
        v = model.get_symbol_by_ontology_term(shared.OXMETA, 'membrane_voltage')
        assert v.type == 'state'

        # Now clamp it to -80mV
        rhs = model.add_number(-80, str(v.units))
        model.set_equation(v, rhs)

        # Check that V is no longer a state
        v = model.get_symbol_by_ontology_term(shared.OXMETA, 'membrane_voltage')
        assert v.type != 'state'

        # TODO: Get dvdt_unit in a more sensible way
        # See: https://github.com/ModellingWebLab/cellmlmanip/issues/133

        # Now make V a state again
        t = model.get_symbol_by_ontology_term(shared.OXMETA, 'time')
        lhs = sp.Derivative(v, t)
        dvdt_units = 'unlikely_unit_name'
        model.add_unit(dvdt_units, [
            {'units': str(v.units)},
            {'units': str(t.units), 'exponent': -1},
        ])
        rhs = model.add_number(0, dvdt_units)
        model.set_equation(lhs, rhs)

        # Check that V is a state again
        v = model.get_symbol_by_ontology_term(shared.OXMETA, 'membrane_voltage')
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
