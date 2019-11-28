import os

import pytest
import sympy as sp

from cellmlmanip import parser
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
    model = load_model_for_session('test_simple_odes')
    return model


@pytest.fixture(scope="session")
def bad_annotation_model():
    model = load_model_for_session('test_bad_annotations')
    return model


