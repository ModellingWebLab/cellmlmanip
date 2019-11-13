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



def test_get_state_symbols():
    model = load_model('aslanidi_model_2009')
    model.get_equation_graph()
    state_symbols = model.get_state_symbols()
    assert len(state_symbols) == 29

