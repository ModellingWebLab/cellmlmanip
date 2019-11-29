import pytest

from cellmlmanip import units


class TestModelUnits:
    def test_symbols(self, model_simple_units):
        symbol = model_simple_units.get_symbol_by_cmeta_id("a")
        assert symbol.is_Symbol
        symbol = model_simple_units.get_symbol_by_cmeta_id("b")
        assert symbol.is_Symbol

    def test_equations(self, model_simple_units):
        symbol_a = model_simple_units.get_symbol_by_cmeta_id("a")
        equation = model_simple_units.get_equations_for([symbol_a])
        assert len(equation) == 1
        assert equation[0].lhs == symbol_a
        assert equation[0].rhs == 1.0

    def test_equations_2(self, model_simple_units):
        symbol_a = model_simple_units.get_symbol_by_cmeta_id("a")
        symbol_b = model_simple_units.get_symbol_by_cmeta_id("b")
        equation = model_simple_units.get_equations_for([symbol_b])
        assert len(equation) == 2
        assert equation[0].lhs == symbol_a
        assert equation[0].rhs == 1.0
        assert equation[1].lhs == symbol_b
        assert equation[1].rhs == 2.0 / symbol_a

    def test_units(self, model_simple_units):
        symbol_a = model_simple_units.get_symbol_by_cmeta_id("a")
        equation = model_simple_units.get_equations_for([symbol_a], strip_units=False)
        assert model_simple_units.units.summarise_units(equation[0].lhs) == 'ms'
        assert model_simple_units.units.summarise_units(equation[0].rhs) == 'ms'

        symbol_b = model_simple_units.get_symbol_by_cmeta_id("b")
        equation = model_simple_units.get_equations_for([symbol_b])
        assert model_simple_units.units.summarise_units(equation[1].lhs) == 'per_ms'
        assert model_simple_units.units.summarise_units(equation[1].rhs) == '1 / ms'
        assert model_simple_units.units.is_unit_equal(model_simple_units.units.summarise_units(equation[1].lhs),
                                                      model_simple_units.units.summarise_units(equation[1].rhs))

    def test_bad_units(self, model_bad_units):
        symbol_a = model_bad_units.get_symbol_by_cmeta_id("a")
        symbol_b = model_bad_units.get_symbol_by_cmeta_id("b")
        equation = model_bad_units.get_equations_for([symbol_b], strip_units=False)
        assert len(equation) == 2
        assert equation[0].lhs == symbol_a
        assert model_bad_units.units.summarise_units(equation[0].lhs) == 'ms'
        with pytest.raises(units.UnitError):
            # cellml file states a (ms) = 1 (ms) + 1 (second)
            model_bad_units.units.summarise_units(equation[0].rhs)

        assert equation[1].lhs == symbol_b
        with pytest.raises(units.UnitError):
            # cellml file states b (per_ms) = power(a (ms), 1 (second))
            model_bad_units.units.summarise_units(equation[1].rhs)
