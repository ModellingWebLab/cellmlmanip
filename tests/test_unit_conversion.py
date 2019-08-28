import pytest
import os
import cellmlmanip

OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata"


@pytest.fixture
def model():
    return cellmlmanip.load_model(os.path.join(os.path.dirname(__file__), 'cellml_files', "test_simple_odes.cellml"))


def test_add_preferred_custom_unit_name(model):
    time_var = model.get_symbol_by_ontology_term(OXMETA, "time")
    assert str(model.units.summarise_units(time_var)) == "ms"
    model.units.add_preferred_custom_unit_name('millisecond', [{'prefix': 'milli', 'units': 'second'}])
    assert str(model.units.summarise_units(time_var)) == "millisecond"
