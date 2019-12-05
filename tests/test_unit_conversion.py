import pytest

from cellmlmanip import units
from . import shared


class TestUnitConversion:
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
