import pytest

from cellmlmanip import units
from . import shared


class TestUnitConversion:
    ###############################################################
    # fixtures

    @pytest.fixture
    def local_model(scope='function'):
        """ Fixture to load a local copy of  the basic_ode model that may get modified. """
        return shared.load_model('basic_ode')

    def test_add_preferred_custom_unit_name(self, simple_ode_model):
        """ Tests Units.add_preferred_custom_unit_name() function. """
        time_var = simple_ode_model.get_symbol_by_ontology_term(shared.OXMETA, "time")
        assert str(simple_ode_model.units.summarise_units(time_var)) == "ms"
        simple_ode_model.units.add_preferred_custom_unit_name('millisecond', [{'prefix': 'milli', 'units': 'second'}])
        assert str(simple_ode_model.units.summarise_units(time_var)) == "millisecond"
        # add_custom_unit does not allow adding already existing units but add_preferred_custom_unit_name does since we
        # cannot know in advance if a model will already have the unit named this way. To test this we add the same unit
        # again
        simple_ode_model.units.add_preferred_custom_unit_name('millisecond', [{'prefix': 'milli', 'units': 'second'}])
        assert str(simple_ode_model.units.summarise_units(time_var)) == "millisecond"

    def test_conversion_factor_original(self, simple_units_model):
        """ Tests Units.get_conversion_factor() function. """
        symbol_b1 = simple_units_model.get_symbol_by_cmeta_id("b_1")
        equation = simple_units_model.get_equations_for([symbol_b1])
        factor = simple_units_model.units.get_conversion_factor(
            quantity=1 * simple_units_model.units.summarise_units(equation[0].lhs),
            to_unit=simple_units_model.units.ureg('us').units)
        assert factor == 1000

    def test_conversion_factor_bad_types(self, simple_units_model):
        """ Tests Units.get_conversion_factor() function for
         cases when arguments are missing or incorrectly typed."""
        symbol_b1 = simple_units_model.get_symbol_by_cmeta_id("b_1")
        equation = simple_units_model.get_equations_for([symbol_b1])
        expression = equation[0].lhs
        to_unit = simple_units_model.units.ureg('us').units
        from_unit = simple_units_model.units.summarise_units(expression)
        quantity = 1 * from_unit
        # no source unit
        with pytest.raises(AssertionError, match='^No unit given as source.*'):
            simple_units_model.units.get_conversion_factor(to_unit=to_unit)
        with pytest.raises(AssertionError, match='^No unit given as source.*'):
            simple_units_model.units.get_conversion_factor(to_unit)

        # no target unit
        with pytest.raises(TypeError):
            simple_units_model.units.get_conversion_factor(from_unit=from_unit)
        # multiple sources
        with pytest.raises(AssertionError, match='^Multiple target.*'):
            simple_units_model.units.get_conversion_factor(to_unit, from_unit=from_unit, quantity=quantity)
        # incorrect types
        with pytest.raises(AssertionError, match='^from_unit must be of type pint:Unit$'):
            simple_units_model.units.get_conversion_factor(to_unit, from_unit=quantity)
        with pytest.raises(AssertionError, match='^quantity must be of type pint:Quantity$'):
            simple_units_model.units.get_conversion_factor(to_unit, quantity=from_unit)
        with pytest.raises(AssertionError, match='^expression must be of type Sympy expression$'):
            simple_units_model.units.get_conversion_factor(to_unit, expression=quantity)

        # unit to unit
        assert simple_units_model.units.get_conversion_factor(to_unit=to_unit, from_unit=from_unit) == 1000
        # quantity to unit
        assert simple_units_model.units.get_conversion_factor(to_unit=to_unit, quantity=quantity) == 1000
        # expression to unit
        assert simple_units_model.units.get_conversion_factor(to_unit=to_unit, expression=expression) == 1000

    def test_conversion_factor_same_units(self, simple_units_model):
        """ Tests Units.get_conversion_factor() function when units are same
        and conversion factor should be '1'. """
        symbol_b = simple_units_model.get_symbol_by_cmeta_id("b")
        equation = simple_units_model.get_equations_for([symbol_b])
        expression = equation[1].rhs
        to_unit = simple_units_model.units.ureg('per_ms').units
        from_unit = simple_units_model.units.summarise_units(expression)
        quantity = 1 * from_unit
        # quantity to unit
        assert simple_units_model.units.get_conversion_factor(to_unit=to_unit, quantity=quantity) == 1
        # unit to unit
        assert simple_units_model.units.get_conversion_factor(to_unit=to_unit, from_unit=from_unit) == 1
        # expression to unit
        assert simple_units_model.units.get_conversion_factor(to_unit=to_unit, expression=expression) == 1

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

    def test_add_input_invalid_args(self, local_model):
        mV_unit = local_model.get_units('mV')
        # name does not exist in model
        with pytest.raises(KeyError):
            local_model.add_input('nonsense_name', mV_unit)

    def test_add_input_state_variable(self, local_model):
        """ Tests the Model.add_input function that changes units. """
        # original state
        def test_original_state(local_model):
            assert len(local_model.variables()) == 3
            symbol_a = local_model.get_symbol_by_cmeta_id('sv11')
            symbol_t = local_model.get_symbol_by_cmeta_id('time')
            assert local_model.get_initial_value(symbol_a) == 2.0
            assert symbol_a.units == 'mV'
            assert symbol_t.units == 'ms'
            assert len(local_model.equations) == 1
            assert str(local_model.equations[0]) == 'Eq(Derivative(_env_ode$sv1, _environment$time), _1.0)'
            return True

        mV_unit = local_model.get_units('mV')
        volt_unit = local_model.get_units('volt')

        assert test_original_state(local_model)
        # test no change in units
        local_model.add_input('env_ode$sv1', mV_unit)
        assert test_original_state(local_model)

        # non-existent unit
        # TODO what if unit not in model

        # change mV to V
        local_model.add_input('env_ode$sv1', volt_unit)
        assert len(local_model.variables()) == 5
        symbol_a = local_model.get_symbol_by_cmeta_id('sv11')
        assert local_model.get_initial_value(symbol_a) == 0.002
        assert symbol_a.units == 'volt'
        assert symbol_a.name == 'env_ode$sv1_converted'
        symbol_t = local_model.get_symbol_by_cmeta_id('time')
        assert symbol_t.units == 'ms'
        symbol_orig = local_model.get_symbol_by_name('env_ode$sv1')
        assert symbol_orig.units == 'mV'
        symbol_derv = local_model.get_symbol_by_name('env_ode$sv1_orig_deriv')
        assert symbol_derv.units == 'mV / ms'
        assert len(local_model.equations) == 3
        assert str(local_model.equations[0]) == 'Eq(_env_ode$sv1, 1000.0*_env_ode$sv1_converted)'
        assert str(local_model.equations[1]) == 'Eq(_env_ode$sv1_orig_deriv, _1.0)'
        assert str(local_model.equations[2]) == 'Eq(Derivative(_env_ode$sv1_converted, _environment$time), ' \
                                                '0.001*_env_ode$sv1_orig_deriv)'

    def test_add_input_free_variable(self, local_model):
        """ Tests the Model.add_input function that changes units. """
        # original state
        def test_original_state(local_model):
            assert len(local_model.variables()) == 3
            symbol_a = local_model.get_symbol_by_cmeta_id('sv11')
            symbol_t = local_model.get_symbol_by_cmeta_id('time')
            assert local_model.get_initial_value(symbol_a) == 2.0
            assert symbol_a.units == 'mV'
            assert symbol_t.units == 'ms'
            assert len(local_model.equations) == 1
            assert str(local_model.equations[0]) == 'Eq(Derivative(_env_ode$sv1, _environment$time), _1.0)'
            return True

        ms_unit = local_model.get_units('ms')
        second_unit = local_model.get_units('second')

        assert test_original_state(local_model)
        # test no change in units
        local_model.add_input('env_ode$time', ms_unit)
        assert test_original_state(local_model)

        # change ms to s
        local_model.add_input('environment$time', second_unit)
        assert len(local_model.variables()) == 5
        symbol_a = local_model.get_symbol_by_cmeta_id('sv11')
        assert local_model.get_initial_value(symbol_a) == 2.0
        assert symbol_a.units == 'mV'
        assert symbol_a.name == 'env_ode$sv1'
        symbol_t = local_model.get_symbol_by_cmeta_id('time')
        assert symbol_t.units == 'second'
        assert symbol_t.name == 'environment$time_converted'
        symbol_orig = local_model.get_symbol_by_name('env_ode$sv1')
        assert symbol_orig.units == 'mV'
        symbol_derv = local_model.get_symbol_by_name('env_ode$sv1_orig_deriv')
        assert symbol_derv.units == 'mV / ms'
        assert len(local_model.equations) == 3
        assert str(local_model.equations[0]) == 'Eq(_environment$time, 1000.0*_environment$time_converted)'
        assert str(local_model.equations[1]) == 'Eq(_env_ode$sv1_orig_deriv, _1.0)'
        assert str(local_model.equations[2]) == 'Eq(Derivative(_env_ode$sv1, _environment$time_converted), ' \
                                                '1000.0*_env_ode$sv1_orig_deriv)'
