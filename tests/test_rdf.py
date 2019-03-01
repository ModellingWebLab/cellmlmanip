import os
import pytest
import sympy

import cellmlmanip
import cellmlmanip.rdf


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
    oxmeta = 'https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#'
    var = model.get_symbol_by_ontology_term(oxmeta, 'time')
    assert isinstance(var, sympy.Symbol)
    assert var.name == 'environment$time'

    # Test getting a non-existent variable
    with pytest.raises(KeyError):
        model.get_symbol_by_ontology_term(oxmeta, 'bert')
    with pytest.raises(KeyError):
        model.get_symbol_by_ontology_term('http://example.com#', 'time')

    # Test bad annotations
    model = load_model('test_bad_annotations')

    # Two variables with the same ID
    with pytest.raises(ValueError) as e:
        model.get_symbol_by_ontology_term(oxmeta, 'time')
    assert 'Multiple variables annotated with' in str(e)

    # Annotation with a cmeta id that doesn't exist
    with pytest.raises(KeyError) as e:
        model.get_symbol_by_ontology_term(oxmeta, 'membrane_potential')
    assert 'No variable with cmeta id' in str(e)

    # Annotation of something that isn't a variable
    with pytest.raises(KeyError) as e:
        model.get_symbol_by_ontology_term(
            oxmeta, 'membrane_fast_sodium_current')
    assert 'No variable with cmeta id' in str(e)

    # Non-local annotation
    # TODO: Add support to allow non-local (but valid, i.e. referring to the
    #       current model) references.
    with pytest.raises(NotImplementedError) as e:
        model.get_symbol_by_ontology_term(
            oxmeta, 'membrane_persistent_sodium_current')
    assert 'Non-local annotations' in str(e)
