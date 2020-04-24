# Module with cellmlmanip tests
import os

import cellmlmanip


OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"


def load_model(name, unit_store=None):
    """ Parses and returns one of the CellML test models """
    if name[-7:] != '.cellml':
        name += '.cellml'
    return cellmlmanip.load_model(
        os.path.join(os.path.dirname(__file__), 'cellml_files', name),
        unit_store=unit_store)
