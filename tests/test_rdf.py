import os
import pytest
import sympy

import cellmlmanip
import cellmlmanip.rdf


@pytest.fixture
def model():
    """ Parses and returns test_simple_odes.cellml. """
    return cellmlmanip.load_model(os.path.join(
        os.path.dirname(__file__), "cellml_files", "test_simple_odes.cellml"))


def test_create_rdf_node():
    # Tests rdf.create_rdf_node

    node = cellmlmanip.rdf.create_rdf_node('http://example.com/', 'hi')
    assert str(node) == 'http://example.com/hi'

    node = cellmlmanip.rdf.create_rdf_node('http://example.com#', 'hi')
    assert str(node) == 'http://example.com#hi'

    node = cellmlmanip.rdf.create_rdf_node('http://example.com', 'hi')
    assert str(node) == 'http://example.com#hi'


def test_get_symbol_by_ontology_term(model):
    # Tests model.get_symbol_by_ontology_term

    # Test getting the time variable
    oxmeta = 'https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#'
    var = model.get_symbol_by_ontology_term(oxmeta, 'time')
    assert isinstance(var, sympy.Symbol)
    assert var.name == 'environment$time'

    # Test getting a non-existent variable
    with pytest.raises(KeyError):
        model.get_symbol_by_ontology_term(oxmeta, 'bert')
    with pytest.raises(KeyError):
        model.get_symbol_by_ontology_term('http://example.com#', 'time')
