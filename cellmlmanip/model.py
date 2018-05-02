"""
Classes representing a CellML model and its components
"""
from io import StringIO

import rdflib


class Component(object):
    """Represents a <component> element in CellML <model>"""

    def __init__(self, name):
        self.name = name
        self.variables = {}
        self.equations = []
        self.numbers = {}

    def __str__(self):
        return "%s\n\tequations: %s\n\tvariables: %s\n\tnumbers: %s" % (self.name,
                                                                        self.equations,
                                                                        self.variables,
                                                                        self.numbers)

    def collect_variable_attributes(self, metadata):
        """
        Called after adding equations and variables from the CellML model
        Collects the disparate bits of information about the equation into a single place
        """
        if not self.equations:
            return

        # Get all unique symbols in the equations
        dummy_symbols = set().union(*[e.free_symbols for e in self.equations])

        # For each symbol in these equations
        for symbol in dummy_symbols:
            # if the symbol is one that's defined as a <variable> in the component
            if symbol.name in self.variables:
                self.variables[symbol.name]['sympy.Dummy'] = symbol
                if symbol in metadata:
                    self.variables[symbol.name].update(metadata[symbol])

            # create variable entry for dummified derivative
            if symbol in metadata:
                if 'sympy.Derivative' in metadata[symbol]:
                    derivative = metadata[symbol]['sympy.Derivative']
                    y_symbol = derivative.free_symbols.pop()
                    x_symbol = derivative.variables[0]
                    # the bound and wrt symbols should be <variable>s in the component
                    self.variables[y_symbol.name]['sympy.Dummy'] = y_symbol
                    self.variables[x_symbol.name]['sympy.Dummy'] = x_symbol
                    # add the dummified derivative symbol itself
                    self.variables[symbol.name] = metadata[symbol]
                    self.variables[symbol.name]['name'] = symbol.name
                elif 'sympy.Number' in metadata[symbol]:
                    self.numbers[symbol] = metadata[symbol]


class Model(object):
    """
    Holds all information about a CellML model and exposes it for manipulation (intention!)
    """
    def __init__(self, name):
        self.name = name
        self.units = {}
        self.components = {}
        self.connections = []
        self.rdf = rdflib.Graph()

    def add_unit(self, units_name: str, unit_elements: dict):
        """
        Adds information about <units> in <model>
        """
        self.units[units_name] = unit_elements

    def add_component(self, component: Component):
        """
        Adds name to list of <component>s in the <model>
        """
        self.components[component.name] = component

    def add_rdf(self, rdf: str):
        """
        Takes RDF string and stores it in an RDFlib.Graph for the model. Can be called repeatedly.
        """
        self.rdf.parse(StringIO(rdf))

    def find_variable(self, search_dict):
        """
        Finds variables in all components in this model that match the given key: value pairs
        :param search_dict:
        :return: a dictionary of key: value pairs
        """
        # a list to collect all matched variables
        matches = []

        # for each component in this model
        for component in self.components.values():
            # for each defined variable in the component
            for variable_attr in component.variables.values():
                # Every search (key, value) pair needs to be in the variable's attributes
                matched = True
                for search_key, search_value in search_dict.items():
                    # If the key is not in the variable or the value is different
                    if search_key not in variable_attr or search_value != variable_attr[search_key]:
                        # This variable doesn't match, break
                        matched = False
                        break
                if matched:
                    # All search (key: value)s were found in this variable
                    matches.append(variable_attr)
        return matches
