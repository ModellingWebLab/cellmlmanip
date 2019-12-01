""" Configuration file for tests setting up the models used as fixtures."""
import pytest

from . import shared


@pytest.fixture(scope='session')
def simple_ode_model():
    """ Returns the test_simple_odes.cellml model for use by testing session.
    Note: do not use if the test attempts to modify the model.
    """
    model = shared.load_model('test_simple_odes')
    return model


@pytest.fixture(scope='session')
def bad_annotation_model():
    """ Returns the test_bad_annotations.cellml model for use by testing session.
    Note: do not use if the test attempts to modify the model.
    """
    model = shared.load_model('test_bad_annotations')
    return model


@pytest.fixture(scope='session')
def basic_model():
    """ Returns the basic_ode.cellml model for use by testing session.
    Note: do not use if the test attempts to modify the model.
    """
    model = shared.load_model('basic_ode')
    return model


@pytest.fixture(scope='session')
def aslanidi_model():
    """ Returns the aslanidi_model_2009.cellml model for use by testing session.
    Note: do not use if the test attempts to modify the model.
    """
    model = shared.load_model('aslanidi_model_2009')
    return model


@pytest.fixture(scope='session')
def hh_model():
    """ Returns the hodgkin_huxley_squid_axon_model_1952_modified.cellml model for use by testing session.
    Note: do not use if the test attempts to modify the model.
    """
    model = shared.load_model('hodgkin_huxley_squid_axon_model_1952_modified')
    return model


@pytest.fixture(scope='session')
def simple_units_model():
    """ Returns the simple_model_units.cellml model for use by testing session.
    Note: do not use if the test attempts to modify the model.
    """
    model = shared.load_model('simple_model_units')
    return model


@pytest.fixture(scope='session')
def bad_units_model():
    """ Returns the simple_model_invalid_units.cellml model for use by testing session.
    Note: do not use if the test attempts to modify the model.
    """
    model = shared.load_model('simple_model_invalid_units')
    return model


