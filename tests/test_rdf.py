import os
import pytest

import cellmlmanip
import cellmlmanip.rdf


@pytest.fixture
def model():
    """ Parses and returns test_simple_odes.cellml. """
    return cellmlmanip.load_model(os.path.join(
        os.path.dirname(__file__), "cellml_files", "test_simple_odes.cellml"))


def test_create_rdf_node():
    # Tests rdf.create_rdf_node
    #TODO
    pass


def test_get_symbol_by_ontology_term(model):
    # Tests model.get_symbol_by_ontology_term

    oxmeta = 'https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#'
    var = model.get_symbol_by_ontology_term(oxmeta, 'time')

    assert var == 'hello'

