"""Main entry-points and functions for interacting with the cellmlmanip library"""
from cellmlmanip.parser import Parser


def load_model(path):
    """Loads a cellml model."""
    model = Parser(path).parse()
    model.get_equation_graph()
    return model
