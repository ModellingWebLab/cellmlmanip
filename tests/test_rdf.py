import os

import pytest
import sympy
import rdflib

import cellmlmanip
import cellmlmanip.rdf


OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"
PYCMLMETA = 'https://chaste.comlab.ox.ac.uk/cellml/ns/pycml#'


def load_model(name):
    """ Parses and returns one of the CellML test models """
    if name[-7:] != '.cellml':
        name += '.cellml'
    return cellmlmanip.load_model(os.path.join(
        os.path.dirname(__file__), 'cellml_files', name))


@pytest.fixture(scope="module")
def test_simple_odes():
    model = load_model('test_simple_odes')
    return model


@pytest.fixture(scope="module")
def test_bad_annotations():
    model = load_model('test_bad_annotations')
    return model


def test_create_rdf_node():
    # Tests rdf.create_rdf_node

    node = cellmlmanip.rdf.create_rdf_node(('http://example.com/', 'hi'))
    assert str(node) == 'http://example.com/hi'
    assert node == rdflib.URIRef('http://example.com/hi')

    node = cellmlmanip.rdf.create_rdf_node(('http://example.com#', 'hi'))
    assert str(node) == 'http://example.com#hi'
    assert node == rdflib.URIRef('http://example.com#hi')

    node = cellmlmanip.rdf.create_rdf_node(('http://example.com', 'hi'))
    assert str(node) == 'http://example.com#hi'
    assert node == rdflib.URIRef('http://example.com#hi')

    node = cellmlmanip.rdf.create_rdf_node('#example')
    assert str(node) == '#example'
    assert node == rdflib.URIRef('#example')

    node = cellmlmanip.rdf.create_rdf_node('example')
    assert str(node) == 'example'
    assert node == rdflib.Literal('example')

    node = cellmlmanip.rdf.create_rdf_node(None)
    assert node is None

    node = rdflib.URIRef('http://example.com#hi')
    assert cellmlmanip.rdf.create_rdf_node(node) == node

    node = rdflib.Literal('example')
    assert cellmlmanip.rdf.create_rdf_node(node) == node


def test_get_symbol_by_ontology_term(test_simple_odes, test_bad_annotations):
    # Tests model.get_symbol_by_ontology_term

    # Test getting the time variable
    model = test_simple_odes
    var = model.get_symbol_by_ontology_term(OXMETA, 'time')
    assert isinstance(var, sympy.Symbol)
    assert var.name == 'environment$time'

    # Test getting a non-existent variable
    with pytest.raises(KeyError):
        model.get_symbol_by_ontology_term(OXMETA, 'bert')
    with pytest.raises(KeyError):
        model.get_symbol_by_ontology_term('http://example.com#', 'time')

    # Test bad annotations
    model = test_bad_annotations

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


def test_get_ontology_terms_by_symbol(test_bad_annotations):
    # Test bad annotations
    model = test_bad_annotations

    # Get v3 from the model, as it does not have cmeta_id, to test this part of the code
    equation_graph = model.get_equation_graph()
    for variable in equation_graph:
        if str(variable) == '_c$v3':
            annotations = model.get_ontology_terms_by_symbol(variable, OXMETA)
            assert len(annotations) == 0


def test_has_ontology_term_by_symbol(test_bad_annotations):
    # Test bad annotations
    model = test_bad_annotations

    # Get v3 from the model, as it does not have cmeta_id, to test this part of the code
    equation_graph = model.get_equation_graph()
    for variable in equation_graph:
        if str(variable) == '_c$v3':
            assert not model.has_ontology_annotation(variable, OXMETA)


def test_get_rdf_annotation(test_simple_odes):
    # Test rdf annotation
    model = test_simple_odes

    named_attributes = []
    named_attrs = model.get_rdf_annotations(subject=model.rdf_identity, predicate=(PYCMLMETA, 'named-attribute'))
    for s, p, attr in named_attrs:
        name = model.get_rdf_value(subject=attr, predicate=(PYCMLMETA, 'name'))
        value = model.get_rdf_value(subject=attr, predicate=(PYCMLMETA, 'value'))
        named_attributes.append({'name': name, 'value': value})

    assert len(named_attributes) == 1
    assert named_attributes[0]['name'] == 'SuggestedForwardEulerTimestep'
    assert named_attributes[0]['value'] == '0.0002'

    # subject as string
    named_attrib = []
    named_attrs = model.get_rdf_annotations(subject='#test_simple_odes', predicate=(PYCMLMETA, 'named-attribute'))
    for s, p, attr in named_attrs:
        name = model.get_rdf_value(subject=attr, predicate=(PYCMLMETA, 'name'))
        value = model.get_rdf_value(subject=attr, predicate=(PYCMLMETA, 'value'))
        named_attrib.append({'name': name, 'value': value})

    assert named_attributes == named_attrib

    params = model.get_symbols_by_rdf((PYCMLMETA, 'modifiable-parameter'), 'yes')
    assert str(params) == '[_single_ode_rhs_const_var$a]'
