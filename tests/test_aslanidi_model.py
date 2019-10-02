import os

import pytest
import sympy

from cellmlmanip import load_model

OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"


class TestAslanidiModel:
    @pytest.fixture(scope="class")
    def model(self):
        cellml = os.path.join(
            os.path.dirname(__file__),
            "cellml_files",
            "aslanidi_model_2009.cellml"
        )
        return load_model(cellml)

    def test_initial_value(self, model):
        graph = model.get_equation_graph()
        membrane_capacitance = model.get_symbol_by_ontology_term(OXMETA, "membrane_capacitance")
        with pytest.raises(KeyError):
            initial_value = model.get_initial_value(membrane_capacitance)
