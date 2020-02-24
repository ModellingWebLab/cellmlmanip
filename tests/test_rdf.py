import pytest
import rdflib
import sympy

import cellmlmanip
import cellmlmanip.rdf
from . import shared

PYCMLMETA = 'https://chaste.comlab.ox.ac.uk/cellml/ns/pycml#'


@pytest.fixture
def local_model(scope='module'):
    """ Local test fixture used for adding rdf. """
    return shared.load_model('test_simple_odes')


def test_create_rdf_node():
    """ Tests rdf.create_rdf_node() function. """

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


def test_get_symbol_by_ontology_term(simple_ode_model, bad_annotation_model):
    """ Tests Model.get_symbol_by_ontology_term() function. """
    # Test getting the time variable
    model = simple_ode_model
    var = model.get_symbol_by_ontology_term((shared.OXMETA, 'time'))
    assert isinstance(var, sympy.Symbol)
    assert var.name == 'environment$time'

    # Test getting a non-existent variable
    with pytest.raises(KeyError):
        model.get_symbol_by_ontology_term((shared.OXMETA, 'bert'))
    with pytest.raises(KeyError):
        model.get_symbol_by_ontology_term(('http://example.com#', 'time'))

    # Test bad annotations
    model = bad_annotation_model

    # Two variables with the same ID
    with pytest.raises(ValueError, match='Multiple variables annotated with'):
        model.get_symbol_by_ontology_term((shared.OXMETA, 'time'))

    # Annotation with a cmeta id that doesn't exist
    with pytest.raises(KeyError, match='No variable with cmeta id'):
        model.get_symbol_by_ontology_term((shared.OXMETA, 'membrane_potential'))

    # Annotation of something that isn't a variable
    with pytest.raises(KeyError, match='No variable with cmeta id'):
        model.get_symbol_by_ontology_term(
            (shared.OXMETA, 'membrane_fast_sodium_current'))

    # Non-local annotation
    # TODO: Add support to allow non-local (but valid, i.e. referring to the
    #       current model) references.
    with pytest.raises(NotImplementedError, match='Non-local annotations'):
        model.get_symbol_by_ontology_term(
            (shared.OXMETA, 'membrane_persistent_sodium_current'))


def test_get_ontology_terms_by_symbol(bad_annotation_model):
    """ Tests Model.get_symbol_by_ontology_term() function when the annotation is not correct. """
    # Test bad annotations
    model = bad_annotation_model

    # Get v3 from the model, as it does not have cmeta_id, to test this part of the code
    for variable in model.graph:
        if str(variable) == '_c$v3':
            annotations = model.get_ontology_terms_by_symbol(variable, shared.OXMETA)
            assert len(annotations) == 0


def test_get_ontology_terms_by_symbol2(hh_model):
    """ Tests Model.get_symbol_by_ontology_term() function when the annotation is not correct. """
    membrane_voltage_var = hh_model.get_symbol_by_ontology_term((shared.OXMETA, "membrane_voltage"))
    annotation = hh_model.get_ontology_terms_by_symbol(membrane_voltage_var, shared.OXMETA)
    assert len(annotation) == 1
    assert annotation[0] == "membrane_voltage"

    # Repeat test with wrong namespace
    annotation = hh_model.get_ontology_terms_by_symbol(membrane_voltage_var, 'http://www.nottingam.ac.uk/')
    assert len(annotation) == 0

    # Repeat test without specifying namespace
    annotation = hh_model.get_ontology_terms_by_symbol(membrane_voltage_var)
    assert len(annotation) == 1
    assert annotation[0] == "membrane_voltage"

    # Get a symbol without cmeta_id and check that get_ontology_terms_by_symbol works
    symbol = [s for s in hh_model.get_state_symbols() if not s.rdf_identity][0]
    annotation = hh_model.get_ontology_terms_by_symbol(symbol)
    assert len(annotation) == 0


def test_has_cmeta_id(basic_model):
    """Tests Model.has_cmeta_id(). """

    # Test with variable ids
    assert basic_model.has_cmeta_id('time')
    assert not basic_model.has_cmeta_id('pepper')

    # Test with model id
    assert basic_model.has_cmeta_id('test_basic_ode')


def test_has_ontology_term_by_symbol(bad_annotation_model):
    """ Tests Model.has_ontology_annotation() function when the annotation is not correct. """
    # Test bad annotations
    model = bad_annotation_model

    # Get v3 from the model, as it does not have cmeta_id, to test this part of the code
    for variable in model.graph:
        if str(variable) == '_c$v3':
            assert not model.has_ontology_annotation(variable, shared.OXMETA)


def test_has_ontology_annotation(hh_model):
    """ Tests Model.has_ontology_annotation() function. """
    membrane_voltage_var = hh_model.get_symbol_by_ontology_term((shared.OXMETA, "membrane_voltage"))
    assert hh_model.has_ontology_annotation(membrane_voltage_var, shared.OXMETA)

    # Repeat test with wrong namespace
    assert not hh_model.has_ontology_annotation(membrane_voltage_var, 'http://www.nottingam.ac.uk/')

    # Repeat test without specifying namespace
    assert hh_model.has_ontology_annotation(membrane_voltage_var)


def test_get_rdf_annotation(simple_ode_model):
    """ Tests Model.get_rdf_annotations() function. """
    model = simple_ode_model

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


def test_add_rdf(local_model):
    """ Tests the Model.add_rdf() function. """

    model = local_model
    named_attributes = []
    named_attrs = model.get_rdf_annotations(subject=model.rdf_identity, predicate=(PYCMLMETA, 'named-attribute'))
    for s, p, attr in named_attrs:
        name = model.get_rdf_value(subject=attr, predicate=(PYCMLMETA, 'name'))
        value = model.get_rdf_value(subject=attr, predicate=(PYCMLMETA, 'value'))
        named_attributes.append({'name': name, 'value': value})

    assert len(named_attributes) == 1
    assert named_attributes[0]['name'] == 'SuggestedForwardEulerTimestep'
    assert named_attributes[0]['value'] == '0.0002'

    model.add_rdf('<rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" >'
                  '<rdf:Description rdf:about="#test_simple_odes">'
                  '<named-attribute xmlns=\"https://chaste.comlab.ox.ac.uk/cellml/ns/pycml#\">'
                  '<rdf:Description>'
                  '<name rdf:datatype="http://www.w3.org/2000/10/XMLSchema#string">AddedAttribute</name>'
                  '<value rdf:datatype="http://www.w3.org/2000/10/XMLSchema#double">2.5</value>'
                  '</rdf:Description>'
                  '</named-attribute>'
                  '</rdf:Description>'
                  '</rdf:RDF>')

    named_attributes = []
    named_attrs = model.get_rdf_annotations(subject=model.rdf_identity, predicate=(PYCMLMETA, 'named-attribute'))
    for s, p, attr in named_attrs:
        name = model.get_rdf_value(subject=attr, predicate=(PYCMLMETA, 'name'))
        value = model.get_rdf_value(subject=attr, predicate=(PYCMLMETA, 'value'))
        named_attributes.append({'name': name, 'value': value})

    assert len(named_attributes) == 2
    # insertion happens randomly so need to check that the new attribute is listed in either position 0/1
    assert named_attributes[0]['name'] == 'AddedAttribute' or named_attributes[1]['name'] == 'AddedAttribute'
    assert named_attributes[0]['value'] == '2.5' or named_attributes[1]['value'] == '2.5'


def test_add_cmeta_id():
    """ Tests Model.add_cmeta_id(). """

    model = shared.load_model('test_simple_odes')

    # Test on variable without a cmeta id
    v = model.get_symbol_by_name('circle_x$x')
    assert v.rdf_identity is None
    model.add_cmeta_id(v)
    assert v.rdf_identity is not None

    # Test on variable with a cmeta id
    v = model.get_symbol_by_ontology_term((shared.OXMETA, 'time'))
    rid = v.rdf_identity
    assert rid is not None
    model.add_cmeta_id(v)
    assert v.rdf_identity == rid

    # Test on variable with a name that's already in use as a cmeta id
    assert model.has_cmeta_id('time')
    v = model.add_variable('time', 'second')
    model.add_cmeta_id(v)


def test_transfer_cmeta_id():
    """ Tests Model.transfer_cmeta_id(). """

    model = shared.load_model('test_simple_odes')

    # Legal move: from variable with a cmeta id to one without
    v1 = model.get_symbol_by_ontology_term((shared.OXMETA, 'time'))
    v2 = model.get_symbol_by_name('circle_x$x')
    assert v1.rdf_identity is not None
    assert v2.rdf_identity is None
    model.transfer_cmeta_id(v1, v2)
    assert v1.rdf_identity is None
    assert v2.rdf_identity is not None

    # Illegal move: from a variable without a cmeta id
    v3 = model.get_symbol_by_name('circle_x$y')
    assert v1.rdf_identity is None
    assert v3.rdf_identity is None
    with pytest.raises(ValueError):
        model.transfer_cmeta_id(v1, v3)

    # Illegal move: to a variable with a cmeta id
    v4 = model.get_symbol_by_ontology_term((shared.OXMETA, 'cytosolic_calcium_concentration'))
    assert v2.rdf_identity is not None
    assert v4.rdf_identity is not None
    with pytest.raises(ValueError):
        model.transfer_cmeta_id(v2, v4)

