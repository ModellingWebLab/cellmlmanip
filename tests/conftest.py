import os

import pytest

from cellmlmanip import load_model


OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"


def load_model_for_session(name):
    """ Parses and returns one of the CellML test models """
    if name[-7:] != '.cellml':
        name += '.cellml'
    return load_model(os.path.join(
        os.path.dirname(__file__), 'cellml_files', name))


@pytest.fixture(scope="session")
def simple_ode_model():
    """ Returns the test_simple_odes.cellml model for use by testing session.
    Note: do not use if the test attempts to modify the model.
    """
    model = load_model_for_session('test_simple_odes')
    return model


@pytest.fixture(scope="session")
def bad_annotation_model():
    """ Returns the test_bad_annotations.cellml model for use by testing session.
    Note: do not use if the test attempts to modify the model.
    """
    model = load_model_for_session('test_bad_annotations')
    return model


@pytest.fixture(scope="session")
def basic_model():
    """ Returns the basic_ode.cellml model for use by testing session.
    Note: do not use if the test attempts to modify the model.
    """
    model = load_model_for_session('basic_ode')
    return model


@pytest.fixture(scope="session")
def aslanidi_model():
    """ Returns the aslanidi_model_2009.cellml model for use by testing session.
    Note: do not use if the test attempts to modify the model.
    """
    model = load_model_for_session('aslanidi_model_2009')
    return model

