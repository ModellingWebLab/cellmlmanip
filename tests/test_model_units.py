import os

import pytest

from cellmlmanip import parser


class TestModelUnits:
    @pytest.fixture(scope="class")
    def parser_instance(self):
        """Parses example CellML and returns model"""
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "simple_model_units.cellml"
        )
        p = parser.Parser(example_cellml)
        return p

    @pytest.fixture(scope="class")
    def model(self, parser_instance):
        model = parser_instance.parse()
        return model

    def test_symbols(self, model):
        symbol = model.get_symbol_by_cmeta_id("a")
        assert symbol.is_Symbol
        symbol = model.get_symbol_by_cmeta_id("b")
        assert symbol.is_Symbol

    def test_equations(self, model):
        symbol_a = model.get_symbol_by_cmeta_id("a")
        equation = model.get_equations_for([symbol_a])
        assert len(equation) == 1
        assert equation[0].lhs == symbol_a
        assert equation[0].rhs == 1.0

    def test_equations_2(self, model):
        symbol_a = model.get_symbol_by_cmeta_id("a")
        symbol_b = model.get_symbol_by_cmeta_id("b")
        equation = model.get_equations_for([symbol_b])
        assert len(equation) == 2
        assert equation[0].lhs == symbol_a
        assert equation[0].rhs == 1.0
        assert equation[1].lhs == symbol_b
        assert equation[1].rhs == 2.0 / symbol_a

    def test_units(self, model):
        symbol_a = model.get_symbol_by_cmeta_id("a")
        equation = model.get_equations_for([symbol_a], strip_units=False)
        assert model.units.summarise_units(equation[0].lhs) == 'ms'
        assert model.units.summarise_units(equation[0].rhs) == 'ms'

        symbol_b = model.get_symbol_by_cmeta_id("b")
        equation = model.get_equations_for([symbol_b])
        assert model.units.summarise_units(equation[1].lhs) == 'per_ms'
        assert model.units.summarise_units(equation[1].rhs) == '1 / ms'
        assert model.units.is_unit_equal(model.units.summarise_units(equation[1].lhs),
                                         model.units.summarise_units(equation[1].rhs))
