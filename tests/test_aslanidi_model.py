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
        assert(model.get_initial_value(membrane_capacitance) == 0.00005)

    def test_initial_value_voltage(self, model):
        model.get_equation_graph()
        membrane_voltage = model.get_symbol_by_ontology_term(OXMETA, "membrane_voltage")
        assert(model.get_initial_value(membrane_voltage) == -80)

    def test_get_derivative_symbols(self, model):
        derivs = model.get_state_symbols()
        derivs_order = model.get_state_symbols(order_by_order_added=True)
        assert len(derivs) == len(derivs_order) == 29

        sorted_derivs = sorted(derivs, key=lambda deriv: str(deriv))
        sorted_derivs_order = sorted(derivs_order, key=lambda deriv: str(deriv))
        for i in range(len(sorted_derivs)):
            assert sorted_derivs[i] == sorted_derivs_order[i]

        assert str(derivs[16]) == '_intracellular_ion_concentrations$K_i'
        assert str(derivs[17]) == '_intracellular_ion_concentrations$Ca_i'
        assert str(derivs_order[16]) == '_intracellular_ion_concentrations$Ca_i'
        assert str(derivs_order[17]) == '_intracellular_ion_concentrations$K_i'
