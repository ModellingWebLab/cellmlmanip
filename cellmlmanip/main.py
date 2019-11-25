"""Main entry-points and functions for interacting with the cellmlmanip library."""
from cellmlmanip.parser import Parser


def load_model(path):
    """Parses a CellML file and returns a :class:`cellmlmanip.Model`."""
    return Parser(path).parse()

