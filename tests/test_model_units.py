import pytest
import os

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

    def test_symbols(self, parser_instance, model):
        model.get_equation_graph(True)  # set up the graph - it is not automatic
        symbol = model.get_symbol_by_cmeta_id("a")
        assert symbol.is_Symbol
        symbol = model.get_symbol_by_cmeta_id("b")
        assert symbol.is_Symbol

    def test_equations(self, parser_instance, model):
        model.get_equation_graph(True)  # set up the graph - it is not automatic
        symbol_a = model.get_symbol_by_cmeta_id("a")
        equation = model.get_equations_for([symbol_a])
        assert len(equation) == 1
        assert equation[0].lhs == symbol_a
        assert equation[0].rhs == 1.0

    def test_equations_2(self, parser_instance, model):
        model.get_equation_graph(True)  # set up the graph - it is not automatic
        symbol_a = model.get_symbol_by_cmeta_id("a")
        symbol_b = model.get_symbol_by_cmeta_id("b")
        equation = model.get_equations_for([symbol_b])
        assert len(equation) == 2
        assert equation[0].lhs == symbol_a
        assert equation[0].rhs == 1.0
        assert equation[1].lhs == symbol_b
        assert equation[1].rhs == 2.0 / symbol_a

    def test_units(self, parser_instance, model):
        model.get_equation_graph(True)  # set up the graph - it is not automatic
        symbol_a = model.get_symbol_by_cmeta_id("a")
        equation = model.get_equations_for([symbol_a])
        assert model.units.summarise_units(equation[0].lhs) == 'ms'
        # this fails but it should not
        # the number 1 has clearly marked units of 'ms'
        # NOTE it fails in the function UnitCalculator::traverse
        assert model.units.summarise_units(equation[0].rhs) == 'ms'

        symbol_b = model.get_symbol_by_cmeta_id("b")
        equation = model.get_equations_for([symbol_b])
        assert model.units.summarise_units(equation[1].lhs) == 'per_ms'
        # this fails in the function UnitCalculator::traverse
#        assert model.units.summarise_units(equation[1].rhs) == 'per_ms'
