import os

import pytest

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

    def test_initial_value_capacitance(self, model):
        model.get_equation_graph()
        membrane_capacitance = model.get_symbol_by_ontology_term(OXMETA, "membrane_capacitance")
        # was raising KeyError but should not
#        with pytest.raises(KeyError):
#            model.get_initial_value(membrane_capacitance)
        assert(model.get_initial_value(membrane_capacitance) == 0.00005)

    def test_initial_value_voltage(self, model):
        model.get_equation_graph()
        membrane_voltage = model.get_symbol_by_ontology_term(OXMETA, "membrane_voltage")
        assert(model.get_initial_value(membrane_voltage) == -80)
