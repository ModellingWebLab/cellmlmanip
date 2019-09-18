import os

import pytest
import sympy

from cellmlmanip import load_model

OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"


class TestSimpleModelUnitsl:
    @pytest.fixture(scope="class")
    def model(self):
        simple_model_units_cellml = os.path.join(
            os.path.dirname(__file__),
            "cellml_files",
            "simple_model_units.cellml"
        )
        return load_model(simple_model_units_cellml)

    def test_get_equations_for(self, model):
        symbol_a = model.get_symbol_by_cmeta_id("a")
        symbol_b = model.get_symbol_by_cmeta_id("b")
        equations = model.get_equations_for([symbol_a, symbol_b], sort_by_input_symbols=False)
        ordered_equation = model.get_equations_for([symbol_a, symbol_b], sort_by_input_symbols=True)
        assert(len(ordered_equation) == len(ordered_equation) == 2)

        ref_eq = [sympy.Eq(sympy.Dummy('A$a'), sympy.numbers.Float(1.0)),
                  sympy.Eq(sympy.Dummy('A$b'), sympy.numbers.Float(2.0) / sympy.Dummy('A$a'))]
        ordered_ref_eq = ref_eq

        # Check equations against expected equations
        for i in range(len(equations)):
            assert str(equations[i]) == str(ref_eq[i])

        # Check equations against expected equations
        for i in range(len(ordered_equation)):
            assert str(ordered_equation[i]) == str(ordered_ref_eq[i])
