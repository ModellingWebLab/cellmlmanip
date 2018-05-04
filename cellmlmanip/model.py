"""
Classes representing a CellML model and its components
"""
import logging
from collections import deque
from io import StringIO

import rdflib
import sympy


class Component(object):
    """Represents a <component> element in CellML <model>"""

    def __init__(self, name):
        self.name = name
        self.variables = {}
        self.equations = []
        self.numbers = {}

    def __str__(self):
        return "Component(\n\tname: %s\n\tequations: %s\n\tvariables: %s\n\tnumbers: %s\n)" % \
               (self.name, self.equations, self.variables, self.numbers)

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

    def make_connections(self):
        """
        Uses public/private interface attributes of variables and model connections to assign
        target variables to their source
        """

        # At this stage we are going to assign variables to new sympy Dummy placeholders
        # This may change in the future once we're confident that everything is working e.g.
        # manipulate the equations directly

        # First, we assign new sympy.Dummy variables for those CellML <variable>s that are only
        # exposed to other components i.e. do not have an "in" value on public/private_interface
        # For each component in the model
        for _, component in self.components.items():
            # For each CellML <variable> in the component
            for _, var_attr in component.variables.items():
                # If this variable does not get its value from another component. There are two
                # possibilities. Either it has:
                # (i) public_interface="out" or
                # (ii) private_interface="out" without public_interface="in"
                if (('public_interface', 'out') in var_attr.items()) or \
                        (('private_interface', 'out') in var_attr.items()
                         and ('public_interface', 'in') not in var_attr.items()):
                    # If it doesn't have a dummy symbol
                    if 'sympy.Dummy' not in var_attr:
                        # This variable was not used in any equations - create a new dummy symbol
                        var_attr['assignment'] = sympy.Dummy(var_attr['name'])
                    else:
                        # The variable is used in an equation & we use the same symbol
                        var_attr['assignment'] = var_attr['sympy.Dummy']

        # Second, we loop over all model connections and create connections between variables

        # Put the connection in a LIFO queue
        connections_to_process = deque(self.connections)

        # For testing: shuffle the order of connections
        # TODO: REMOVE!
        from random import shuffle
        shuffle(connections_to_process)

        # While we still have connections left to process
        while connections_to_process:
            # Get connection at front of queue
            connection = connections_to_process.popleft()
            logging.info("Try to connect %s and %s", *connection)
            success = self.__connect(connection)
            if success:
                logging.info('Connected.')
            else:
                logging.info('Cannot connect (source does not have assignment).')
                connections_to_process.append(connection)
            # TODO: track looping and break if we can't exit

    def __get_connection_parts(self, connection):
        ((component_1, variable_1), (component_2, variable_2)) = connection
        variable_1_attributes = self.components[component_1].variables[variable_1]
        variable_2_attributes = self.components[component_2].variables[variable_2]
        return component_1, variable_1_attributes, component_2, variable_2_attributes

    def __connect(self, connection):
        comp_1_name, var_1_attr, comp_2_name, var_2_attr = self.__get_connection_parts(connection)
        # Determine the source and target variables
        if (('public_interface', 'out') in var_1_attr.items() or
                ('private_interface', 'out') in var_1_attr.items()) and (
                    ('public_interface', 'in') in var_2_attr.items() or
                    ('private_interface', 'in') in var_2_attr.items()):
            return self.__connect_with_direction(comp_1_name, var_1_attr, comp_2_name, var_2_attr)
        elif (('public_interface', 'out') in var_2_attr.items() or
              ('private_interface', 'out') in var_2_attr.items()) and (
                  ('public_interface', 'in') in var_1_attr.items() or
                  ('private_interface', 'in') in var_1_attr.items()):
            return self.__connect_with_direction(comp_2_name, var_2_attr, comp_1_name, var_1_attr)
        raise RuntimeError("Cannot determine the source & target for connection %s" % connection)

    def __connect_with_direction(self, source_component, source_variable,
                                 target_component, target_variable):
        logging.info('    Source: %s -> %s', source_component, source_variable)
        logging.info('    Target: %s -> %s', target_component, target_variable)
        # If the source variable has already been assigned a final symbol
        if 'assignment' in source_variable:
            # If source/target variable is in the same unit
            if source_variable['units'] == target_variable['units']:
                # Direct substitution is possible
                target_variable['assignment'] = source_variable['assignment']
            else:
                # Requires a conversion, so we add an equation to the component that assigns the
                # target dummy variable to the source variable (unit conversion handled separately)
                self.components[target_component].equations.append(
                    sympy.Eq(target_variable['sympy.Dummy'], source_variable['assignment'])
                )
                logging.info('    New target eq: %s -> %s',
                             target_component, self.components[target_component].equations[-1])

                # The assigned symbol for this variable is itself
                target_variable['assignment'] = target_variable['sympy.Dummy']
            logging.info('    Updated target: %s -> %s', target_component, target_variable)
            return True
        return False

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
