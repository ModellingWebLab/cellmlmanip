import pytest
import os

from cellmlmanip import parser


class TestModelUnits():
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
        symbolA = model.get_symbol_by_cmeta_id("a")
        equation = model.get_equations_for([symbolA])
        assert len(equation) == 1
        assert equation[0].lhs == symbolA
        assert equation[0].rhs == 1.0

        symbolB = model.get_symbol_by_cmeta_id("b")
        # Note calling this function a second time adds to the list of equations
        # I'm not sure that is expected behaviour -
        # I have only asked for equations for 'b'
        equation1 = model.get_equations_for([symbolB])
        assert len(equation1) == 2
        assert equation1[1].lhs == symbolB
        assert equation1[1].rhs == 2.0 / symbolA

    def test_units(self, parser_instance, model):
        model.get_equation_graph(True)  # set up the graph - it is not automatic
        symbolA = model.get_symbol_by_cmeta_id("a")
        equation = model.get_equations_for([symbolA])
        assert model.units.summarise_units(equation[0].lhs) == 'ms'
        # this fails but it should not
        # the number 1 has clearly marked units of 'ms'
        # NOTE it fails in the function UnitCalculator::traverse
#        assert model.units.summarise_units(equation[0].rhs) == 'ms'

        symbolB = model.get_symbol_by_cmeta_id("b")
        equation = model.get_equations_for([symbolB])
        assert model.units.summarise_units(equation[1].lhs) == 'per_ms'
        # this fails in the function UnitCalculator::traverse
#        assert model.units.summarise_units(equation[1].rhs) == 'per_ms'
