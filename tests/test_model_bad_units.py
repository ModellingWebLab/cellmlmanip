import os

import pytest

from cellmlmanip import parser, units


class TestModelBadUnits:
    @pytest.fixture(scope="class")
    def parser_instance(self):
        """Parses example CellML and returns model"""
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "simple_model_invalid_units.cellml"
        )
        p = parser.Parser(example_cellml)
        return p

    @pytest.fixture(scope="class")
    def model(self, parser_instance):
        model = parser_instance.parse()
        return model

    def test_units(self, parser_instance, model):
        symbol_a = model.get_symbol_by_cmeta_id("a")
        symbol_b = model.get_symbol_by_cmeta_id("b")
        equation = model.get_equations_for([symbol_b], strip_units=False)
        assert len(equation) == 2
        assert equation[0].lhs == symbol_a
        assert model.units.summarise_units(equation[0].lhs) == 'ms'
        with pytest.raises(units.UnitError):
            # cellml file states a (ms) = 1 (ms) + 1 (second)
            model.units.summarise_units(equation[0].rhs)

        assert equation[1].lhs == symbol_b
        with pytest.raises(units.UnitError):
            # cellml file states b (per_ms) = power(a (ms), 1 (second))
            model.units.summarise_units(equation[1].rhs)
