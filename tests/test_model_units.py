import pytest

from cellmlmanip import units


# TODO some tests here are repeats and may not be necessary
class TestModelUnits:
    def test_symbols(self, simple_units_model):
        """ Tests the Model.get_symbol_by_cmeta_id function."""
        symbol = simple_units_model.get_symbol_by_cmeta_id("a")
        assert symbol.is_Symbol
        symbol = simple_units_model.get_symbol_by_cmeta_id("b")
        assert symbol.is_Symbol

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
