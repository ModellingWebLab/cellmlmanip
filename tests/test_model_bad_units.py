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
        model.get_equation_graph(True)  # set up the graph - it is not automatic
        symbol_a = model.get_symbol_by_cmeta_id("a")
        symbol_b = model.get_symbol_by_cmeta_id("b")
        equation = model.get_equations_for([symbol_b])
        assert len(equation) == 2
        assert equation[0].lhs == symbol_a
        assert equation[0].rhs == 2.0
        assert model.units.summarise_units(equation[0].lhs) == 'ms'
        # cellml file states a (ms) = 1 (ms) + 1 (second)
        # this does not raise and error
        # instead it returns 2 dimensionless
        # with pytest.raises(units.UnitError):
        model.units.summarise_units(equation[0].rhs)

        assert equation[1].lhs == symbol_b
        assert equation[1].rhs == symbol_a**1.0
        # cellml file states b (per_ms) = power(a (ms), 1 (second))
        # this does raise an error
        with pytest.raises(units.UnitError):
            model.units.summarise_units(equation[1].rhs)
