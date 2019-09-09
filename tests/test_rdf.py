import os

import pytest
import sympy

import cellmlmanip
import cellmlmanip.rdf

OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"


def load_model(name):
    """ Parses and returns one of the CellML test models """
    if name[-7:] != '.cellml':
        name += '.cellml'
    return cellmlmanip.load_model(os.path.join(
        os.path.dirname(__file__), 'cellml_files', name))


def test_create_rdf_node():
    # Tests rdf.create_rdf_node

    node = cellmlmanip.rdf.create_rdf_node('http://example.com/', 'hi')
    assert str(node) == 'http://example.com/hi'

    node = cellmlmanip.rdf.create_rdf_node('http://example.com#', 'hi')
    assert str(node) == 'http://example.com#hi'

    node = cellmlmanip.rdf.create_rdf_node('http://example.com', 'hi')
    assert str(node) == 'http://example.com#hi'


def test_get_symbol_by_ontology_term():
    # Tests model.get_symbol_by_ontology_term

    # Test getting the time variable
    model = load_model('test_simple_odes')
    var = model.get_symbol_by_ontology_term(OXMETA, 'time')
    assert isinstance(var, sympy.Symbol)
    assert var.name == 'environment$time'

    # Test getting a non-existent variable
    with pytest.raises(KeyError):
        model.get_symbol_by_ontology_term(OXMETA, 'bert')
    with pytest.raises(KeyError):
        model.get_symbol_by_ontology_term('http://example.com#', 'time')

    # Test bad annotations
    model = load_model('test_bad_annotations')

    # Two variables with the same ID
    with pytest.raises(ValueError, match='Multiple variables annotated with'):
        model.get_symbol_by_ontology_term(OXMETA, 'time')

    # Annotation with a cmeta id that doesn't exist
    with pytest.raises(KeyError, match='No variable with cmeta id'):
        model.get_symbol_by_ontology_term(OXMETA, 'membrane_potential')

    # Annotation of something that isn't a variable
    with pytest.raises(KeyError, match='No variable with cmeta id'):
        model.get_symbol_by_ontology_term(
            OXMETA, 'membrane_fast_sodium_current')

    # Non-local annotation
    # TODO: Add support to allow non-local (but valid, i.e. referring to the
    #       current model) references.
    with pytest.raises(NotImplementedError, match='Non-local annotations'):
        model.get_symbol_by_ontology_term(
            OXMETA, 'membrane_persistent_sodium_current')


def test_get_ontology_term_by_symbol():
    # Test bad annotations
    model = load_model('test_bad_annotations')
    v1 = model.get_symbol_by_cmeta_id('v1')

    # Get v3 from the model, as it does not have cmeta_id, to test this part of the code
    equation_graph = model.get_equation_graph()
    for variable in equation_graph:
        if str(variable) == '_c$v3':
            annotations = model.get_ontology_terms_by_symbol(OXMETA, variable)
            assert len(annotations) == 0
