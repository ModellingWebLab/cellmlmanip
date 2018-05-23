"""
Classes representing a CellML model and its components
"""
import logging
from collections import deque
from io import StringIO

import rdflib
import sympy
import sympy.physics.units as units


class Component(object):
    """Represents a <component> element in CellML <model>"""

    def __init__(self, name):
        self.name = name
        self.variables = {}
        self.equations = []
        self.numbers = {}

    def __str__(self):
        return "Component(\n  name: %s\n  equations: %s\n  variables: %s\n  numbers: %s\n)" % \
               (self.name, self.equations, self.variables, self.numbers)

    def collect_variable_attributes(self, metadata):
        """
        Called after adding equations and variables from the CellML model
        Collects the disparate bits of information about the equation into a single place
        """
        if not self.equations:
            return

        # Get all unique symbols in the equations
        dummy_symbols = set()
        for equation in self.equations:
            for free_symbol in equation.free_symbols:
                dummy_symbols.add(free_symbol)
            # Collect bound variables for derivatives (not returned by .free_symbols)
            for derivative in equation.atoms(sympy.Derivative):
                dummy_symbols.add(derivative.variables[0])  # TODO: variables for degree > 1

        # For each of the symbols in these equations
        for symbol in dummy_symbols:
            # If this symbol is a number
            if symbol in metadata and 'sympy.Number' in metadata[symbol]:
                # save the number metadata (units etc.) and skip variable
                self.numbers[symbol] = metadata[symbol]
                continue

            # if the symbol is one that's defined as a <variable> in the component
            variable_name = symbol.name.split('__')[1]
            if symbol.name.split('__')[1] in self.variables:
                self.variables[variable_name]['sympy.Dummy'] = symbol


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

    def __str__(self):
        return '\n'.join([str(v) for k, v in self.components.items()])

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
                        var_attr['sympy.Dummy'] = var_attr['assignment'] = sympy.Dummy(
                            component.name + '__' + var_attr['name'])
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

        # All connections have been made, now substitute the original dummy with new dummy
        for _, component in self.components.items():
            for _, variable in component.variables.items():
                for index, equation in enumerate(component.equations):
                    # If this variable is used in an equation, it will have been set a sympy.Dummy
                    if 'sympy.Dummy' in variable:
                        # If the variable has been assigned [a new dummy symbol]
                        if 'assignment' in variable:
                            # Replace the original dummy with the assign dummy symbol
                            component.equations[index] = equation.subs(
                                {variable['sympy.Dummy']: variable['assignment']})
                        else:
                            variable['assignment'] = variable['sympy.Dummy']
                    # The variable does not have a sympy.Dummy variable set - why??
                    else:
                        logging.warning('Variable (%s) in component (%s) not assigned a dummy',
                                        variable['name'], component.name)

    def __get_connection_endpoints(self, connection):
        ((component_1, variable_1), (component_2, variable_2)) = connection
        variable_1_attributes = self.components[component_1].variables[variable_1]
        variable_2_attributes = self.components[component_2].variables[variable_2]
        return component_1, variable_1_attributes, component_2, variable_2_attributes

    def __connect(self, connection):
        cmp_name_1, var_att_1, cmp_name_2, var_att_2 = self.__get_connection_endpoints(connection)
        # Determine the source and target variables
        if (('public_interface', 'out') in var_att_1.items() or
                ('private_interface', 'out') in var_att_1.items()) and (
                    ('public_interface', 'in') in var_att_2.items() or
                    ('private_interface', 'in') in var_att_2.items()):
            return self.__connect_with_direction(cmp_name_1, var_att_1, cmp_name_2, var_att_2)
        elif (('public_interface', 'out') in var_att_2.items() or
              ('private_interface', 'out') in var_att_2.items()) and (
                  ('public_interface', 'in') in var_att_1.items() or
                  ('private_interface', 'in') in var_att_1.items()):
            return self.__connect_with_direction(cmp_name_2, var_att_2, cmp_name_1, var_att_1)
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
        # The source variable has not been assigned a symbol, so we can't make this connection
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


class QuantityStore(object):
    """
    Holds sympy.physics.unit.Quantity objects for the model. Can be initialised with <units>
    definitions from the CellML <model>. The get_quantity() methods will find the right
    quantity to return, given the name of the unit (from the internal store, Sympy built-in, or
    constructed anew)
    """

    # Aliases for units to their Sympy equivalents
    aliases = {
        'litre': 'liter',
        'litres': 'liters',
        'metre': 'meter',
        'metres': 'meters'
    }

    def __init__(self, cellml_def=None):
        self.store = {}
        self.cellml_definitions = cellml_def if cellml_def else {}
        self.sympify_context = {}

        from sympy.core.compatibility import exec_
        exec_('from sympy.physics.units import *', self.sympify_context)

    def get_quantity(self, unit_name):
        """
        Given the name of the unit, this will either (i) return a Quantity from the internal store
        as its already been resolved previously (ii) return the Quantity from Sympy if it a built-in
        name or (iii) construct and return a new Quantity object using the <units><unit></units>
        definition in the CellML <model>
        """
        # This name is aliased to some other unit name (e.g. british -> american spelling)
        if unit_name in QuantityStore.aliases:
            unit_name = QuantityStore.aliases[unit_name]

        # If we've already sourced the quantity for this name
        if unit_name in self.store:
            return self.store[unit_name]

        # If this unit name is one of the built-in Sympy quantities
        if hasattr(units, unit_name):
            self.store[unit_name] = getattr(units, unit_name)
            return self.store[unit_name]

        # If this unit name is defined in the CellML model
        if unit_name in self.cellml_definitions:
            # Make a Sympy Quantity object for this unit definition
            quantity = self._make_cellml_quantity(unit_name)
            self.store[unit_name] = quantity
            return self.store[unit_name]

        raise RuntimeError('Cannot find the unit with name (%s)' % unit_name)

    def _sympify(self, string):
        # logging.info('sympy.sympify(%s)', string)
        return sympy.sympify(string, locals=self.sympify_context)

    def _make_cellml_quantity(self, name):
        full_unit_expr = []
        full_dimension = []

        # For each of the <unit> elements for this unit definition
        for unit_element in self.cellml_definitions[name]:
            # Source the <unit units="XXX"> from our store
            unit_as_quantity = self.get_quantity(unit_element['units'])

            # Add this unit to the sympify context if necessary
            if unit_element['units'] not in self.sympify_context:
                self.sympify_context[unit_element['units']] = unit_as_quantity

            # Construct a string representing the expression and dimensions for this <unit>
            expr = str(unit_as_quantity.name)
            dimension = str(unit_as_quantity.args[1].args[0])

            if 'prefix' in unit_element:
                expr = '%s * %s' % (unit_element['prefix'], expr)

            expr = '(%s)' % expr

            if 'exponent' in unit_element:
                expr = '%s**%s' % (expr, unit_element['exponent'])
                # Exponent of the <unit> also affects the dimension
                dimension = '(%s)**%s' % (dimension, unit_element['exponent'])

            expr = '(%s)' % expr

            full_unit_expr.append(expr)
            full_dimension.append('(%s)' % dimension)

        # Have sympy evaluate the collected unit dimension and expression
        full_dimension = self._sympify('*'.join(full_dimension))
        full_unit_expr = self._sympify('*'.join(full_unit_expr))

        # Create and return the Quantity object for this CellML <units>
        quantity = units.Quantity(name, full_dimension, full_unit_expr)
        logging.info('%s=Quantity("%s", %s, %s)',
                     name, name, full_dimension.args[0], full_unit_expr)
        return quantity
