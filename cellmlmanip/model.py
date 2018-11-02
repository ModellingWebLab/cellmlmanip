"""
Classes representing a CellML model and its components
"""
import logging
from collections import deque
from io import StringIO
from typing import Dict, List, Tuple

import rdflib
import sympy
import sympy.physics.units as units


class Component(object):
    """Represents a <component> element in CellML <model>"""

    def __init__(self, name: str, model: 'Model') -> None:
        """Create a new Component object

        :param name: the name of the component, usually from <component name="">
        :param model: a reference to the parent Model
        """
        self.name: str = name
        self.model: 'Model' = model
        self.variables: Dict[str, Dict] = {}
        self.equations: List[sympy.Eq] = []
        self.numbers: Dict[sympy.Dummy, Dict] = {}

    def __str__(self):
        """Pretty-print object"""
        return "Component(\n  name: %s\n  equations: %s\n  variables: %s\n  numbers: %s\n)" % (
            self.name, self.equations, self.variables, self.numbers
        )

    def collect_variable_attributes(self, metadata):
        """Called after adding equations and variables from the CellML model and collects the
        disparate bits of information about the equation into a single place. Associates the
        Sympy symbol for an object with its variable definition if possible.

        :param metadata: returned by the mathml2sympy transpiler, containing information that cannot
            be returned in a sympy expression (e.g. the unit of a number)
        """

        # Skip if this component doesn't have any equations
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
            if variable_name in self.variables:
                self.variables[variable_name]['sympy.Dummy'] = symbol

    def add_units_to_equations(self):
        """Inserts unit attributes associated with symbols into equations.

        WARNING: This method replaces/mutates the existing equations i.e. the original equations
        are lost.
        """
        # For each equation in this component
        for index, equation in enumerate(self.equations):
            # extract the unit information from the equation
            unit_info = self._collect_units(equation)

            # Re-purpose the unit_info dict into a substitution dictionary
            # to replace the "symbol-or-number" with "unit * symbol-or-number"
            for key, value in unit_info.items():
                # If this unit information is about a number
                if key in self.numbers and self.numbers[key]['sympy.Number'] != 0:
                    # Replace the dummy number with the sympy.Number (generated by the parser)
                    unit_info[key] = self.numbers[key]['sympy.Number'] * value
                # Otherwise, straightforward insertion of unit
                else:
                    unit_info[key] = key * value

            # Do the substitutions (WARNING: irreversible!)
            self.equations[index] = equation.subs(unit_info)

    def _collect_units(self, expr):
        """Descends into the given Sympy expression and returns the appropriate units (if any)"""
        logging.debug("_collect_units(%s)", expr)

        expr_unit_info = self._find_unit_info(expr)
        if expr_unit_info:
            expr_unit_info = self.model.units.get_quantity(expr_unit_info)
            unit_info = {expr: expr_unit_info}
            logging.debug('Found unit for "%s" ⟶ "%s"', expr, expr_unit_info)
        else:
            unit_info: Dict = {}

        # Traverse down the arguments of this expression, and collect their unit information
        for args in expr.args:
            unit_info.update(self._collect_units(args))

        # Special handling for Derivatives
        if isinstance(expr, sympy.Derivative):
            unit_info = {expr: unit_info[expr.args[0]] / unit_info[expr.args[1]]}
        return unit_info

    def _find_unit_info(self, expr):
        """Takes an expression (part of an equation) and searches in the component variables
        and numbers to get the associated information about this expression/symbol

        :param expr: a Sympy expression. The interesting ones are leaf nodes of the expression
        """
        variable_with_assigned = self.find_variable({'assignment': expr})
        if variable_with_assigned:
            return variable_with_assigned[0]['units']

        variable_with_dummy = self.find_variable({'sympy.Dummy': expr})
        if variable_with_dummy:
            return variable_with_dummy[0]['units']

        for dummy_number, number_info in self.numbers.items():
            if dummy_number == expr:
                return number_info['cellml:units']

        # Note unit information about this expression in this component. It might be referring to
        # a variable in a different component
        variable_with_assigned = self.model.find_variable({'assignment': expr})
        if variable_with_assigned:
            return variable_with_assigned[0]['units']

        # No information found about this
        return None

    def find_variable(self, search_dict):
        """Finds variables in this model that fulfil the search criteria.

        Each model variable has a number of attributes e.g. name, public_interface, assignment,
        cmeta:id etc. Search criteria are described as key-value pairs in a dictionary. Those
        model variables that fulfil _all_ key-value pairs in the search_dict are returned.

        :param search_dict: a dictionary providing the search criteria
        :return: a list of variables in the model matching the search criteria
        """

        # throw error if empty search criteria given
        if search_dict is None or not search_dict:
            raise RuntimeError("Search criteria cannot be empty.")

        # a list to collect all matched variables
        matches = []

        # for each defined variable in the component
        for variable in self.variables.values():
            # Every search (key, value) pair needs to be in the variable's attributes
            matched = True
            for search_key, search_value in search_dict.items():
                # If the key is not in the variable or the value is different
                if search_key not in variable or search_value != variable[search_key]:
                    # This variable doesn't match, break
                    matched = False
                    break
            if matched:
                # All search (key: value)s were found in this variable
                matches.append(variable)
        return matches


class Model(object):
    """Holds all information about a CellML model and exposes it for manipulation (intention!)
    """
    def __init__(self, name: str) -> None:
        """Create a new instance of Model object
        :param name: the name of the model e.g. from <model name="">
        """
        self.name: str = name
        self.units: 'QuantityStore' = None
        self.components: Dict[str, Component] = {}
        self.connections: List[Tuple] = []
        self.rdf: rdflib.Graph = rdflib.Graph()

    def __str__(self):
        """Pretty-print each of the components in this model"""
        return '\n'.join([str(v) for v in self.components.values()])

    def add_unit(self, units_elements: dict):
        """Adds information about <units> in <model>
        """
        self.units = QuantityStore(units_elements)

    def add_component(self, component: Component):
        """Adds name to list of <component>s in the <model>
        """
        if component.name in self.components:
            raise RuntimeError("%s already exists. Check CellML." % component.name)

        self.components[component.name] = component

    def add_rdf(self, rdf: str):
        """ Takes RDF string and stores it in an RDFlib.Graph for the model. Can be called
        repeatedly.
        """
        self.rdf.parse(StringIO(rdf))

    def _is_free_variable(self, variables: Dict):
        """Checks whether a variable gets its value from another component. There are two
        possibilities. Either it has:
            (i) public_interface="out" or
            (ii) private_interface="out" without public_interface="in"
        """
        return (('public_interface', 'out') in variables.items()) or \
               (('private_interface', 'out') in variables.items() and
                ('public_interface', 'in') not in variables.items())

    def make_connections(self):
        """Uses public/private interface attributes of variables and model connections to assign
        target variables to their source
        """

        # At this stage we are going to assign variables to new sympy Dummy placeholders
        # This may change in the future once we're confident that everything is working e.g.
        # manipulate the equations directly

        # First, we assign new sympy.Dummy variables for those CellML <variable>s that are only
        # exposed to other components i.e. do not have an "in" value on public/private_interface
        # For each component in the model
        for component in self.components.values():
            # For each CellML <variable> in the component
            for var_attr in component.variables.values():
                # If this variable does not get its value from another component
                if self._is_free_variable(var_attr):
                    # If it doesn't have a dummy symbol
                    if 'sympy.Dummy' not in var_attr:
                        # This variable was not used in any equations - create a new dummy symbol
                        var_attr['sympy.Dummy'] = var_attr['assignment'] = sympy.Dummy(
                            component.name + '__' + var_attr['name'])
                    else:
                        # The variable is used in an equation & we use the same symbol
                        var_attr['assignment'] = var_attr['sympy.Dummy']

        # Second, we loop over all model connections and create connections between variables

        # Put the connections in a LIFO queue
        connections_to_process = deque(self.connections)

        # For testing: shuffle the order of connections
        # TODO: REMOVE!
        from random import shuffle
        shuffle(connections_to_process)

        # While we still have connections left to process
        while connections_to_process:
            # Get connection at front of queue
            connection = connections_to_process.popleft()
            logging.debug("Try to connect %s and %s", *connection)
            success = self.__connect(connection)
            if success:
                logging.debug('Connected.')
            else:
                logging.debug('Cannot connect (source does not have assignment).')
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
        """Pull out the variable dict of the component for the two endpoints of the connection

        :param connection: single connection tuple as created by the parser
        """
        ((component_1, variable_1), (component_2, variable_2)) = connection
        variable_1_attributes = self.components[component_1].variables[variable_1]
        variable_2_attributes = self.components[component_2].variables[variable_2]
        return component_1, variable_1_attributes, component_2, variable_2_attributes

    def __connect(self, connection):
        """Takes a CellML connection and attempts to resolve the connect by assigning the target
        variable to the assigned source variable

        :param connection: a single connection tuple, created by the CellML parser
            ((component_1, variable_1), (component_2, variable_2))
        """
        cmp_name_1, var_att_1, cmp_name_2, var_att_2 = self.__get_connection_endpoints(connection)

        # Determine which are the source and target variables, and connect from->to
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
        """Given the source and target component and variable, create a connection by assigning
        the symbol from the source to the target. If units are not the same, it will add an equation
        to the component reflecting the relationship. If a symbol has not been assigned to the
        source variable, then return False.

        :param source_component: name of source component
        :param source_variable: name of source variable
        :param target_component: name of target component
        :param target_variable: name of target variable
        """
        logging.debug('    Source: %s ⟶ %s', source_component, source_variable)
        logging.debug('    Target: %s ⟶ %s', target_component, target_variable)
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
                logging.info('    New target eq: %s ⟶ %s', target_component,
                             self.components[target_component].equations[-1])

                # The assigned symbol for this variable is itself
                target_variable['assignment'] = target_variable['sympy.Dummy']
            logging.debug('    Updated target: %s ⟶ %s', target_component, target_variable)
            return True
        # The source variable has not been assigned a symbol, so we can't make this connection
        return False

    def find_variable(self, search_dict):
        """Finds variables in all components in this model that match the given key: value pairs

        :param search_dict:
        :return: a dictionary of key: value pairs
        """
        # a list to collect all matched variables
        matches = []

        # for each component in this model
        for component in self.components.values():
            matches.extend(component.find_variable(search_dict))
        return matches


class QuantityStore(object):
    """Holds sympy.physics.unit.Quantity objects for the model. Can be initialised with <units>
    definitions from the CellML <model>. The get_quantity() methods will find the right
    quantity to return, given the name of the unit (from the internal store, Sympy built-in, or
    constructed anew)
    """

    # Aliases for units to their Sympy equivalent
    UNIT_ALIASES = {
        'litre': 'liter',
        'metre': 'meter',
    }

    # The full list of supported CellML units
    # Taken from https://www.cellml.org/specifications/cellml_1.1/#sec_units
    # Some are not defined by Sympy, see comments
    CELLML_UNITS = [
        # Base SI units
        'ampere',
        'candela',
        'kelvin',
        'kilogram',
        'meter',
        'metre',  # in self.UNIT_ALIASES
        'mole',
        'second',

        # Derived SI units
        'becquerel',  # see __add_custom_units()
        # 'celsius',  # Not supported in CellML 2.0 nor Sympy
        'coulomb',
        'farad',
        'gray',  # see __add_custom_units()
        'henry',
        'hertz',
        'joule',
        'katal',  # see __add_custom_units()
        'lumen',  # see __add_custom_units()
        'lux',
        'newton',
        'ohm',
        'pascal',
        'radian',
        'siemens',
        'sievert',  # see __add_custom_units()
        'steradian',
        'tesla',
        'volt',
        'watt',
        'weber',

        # Convenience units
        'dimensionless',  # see __add_custom_units()
        'gram',
        'liter',
        'litre',  # in self.UNIT_ALIASES
    ]

    def __init__(self, cellml_def=None):
        """Initialise a QuantityStore instance
        :param cellml_def: a dictionary of <units> definitions from the CellML model. See parser
        for format, essentially: {'name of unit': { [ unit attributes ], [ unit attributes ] } }
        """
        # Initialise the store
        self.store = {}

        # Add required quantities not provided by Sympy
        self.__add_custom_units()

        self.cellml_definitions = cellml_def if cellml_def else {}
        self.sympify_context = {}

        # Setup the context with Sympy unit definitions for sympify
        from sympy.core.compatibility import exec_
        exec_('from sympy.physics.units import *', self.sympify_context)

    def __add_custom_units(self):
        """Adds custom Sympy dimensions and quantities that aren't provided by default but
        required by the CellML specification
        """
        self.store['dimensionless'] = units.Quantity('dimensionless',
                                                     units.Dimension(1),
                                                     sympy.Rational(1, 1))

        self.store['becquerel'] = units.Quantity('becquerel', units.Dimension('activity'), 1, 'Bq')

        # Taken from https://github.com/sympy/sympy/pull/13658
        # gray is a J/kg physical quantity
        self.store['gray'] = units.Quantity('gray',
                                            units.energy/units.mass,
                                            units.meter**2/units.second**2)

        self.store['katal'] = units.Quantity('katal',
                                             units.amount_of_substance/units.time,
                                             units.mol/units.second)

        self.store['lumen'] = units.Quantity('lumen',
                                             units.luminous_intensity * units.steradian.dimension,
                                             units.candela * units.steradian)

        # See https://en.wikipedia.org/wiki/Sievert#Definition for relationship with gray
        # sievert is J/kg biological effect
        self.store['sievert'] = units.Quantity('sievert',
                                               units.energy/units.mass,
                                               units.meter**2/units.second**2)

    def get_quantity(self, unit_name):
        """Given the name of the unit, this will either (i) return a Quantity from the internal
        store as its already been resolved previously (ii) return the Quantity from Sympy if it a
        built-in name or (iii) construct and return a new Quantity object using the
        <units><unit></units> definition in the CellML <model>
        """
        # Correct any aliases (e.g. british -> american spelling)
        unit_name = QuantityStore.UNIT_ALIASES.get(unit_name, unit_name)

        # If we've already sourced the quantity for this name
        if unit_name in self.store:
            return self.store[unit_name]

        # If this unit name is defined in the CellML model
        if unit_name in self.cellml_definitions:
            # Make a Sympy Quantity object for this unit definition
            quantity = self._make_cellml_quantity(unit_name)
            self.store[unit_name] = quantity
            return self.store[unit_name]

        # If this unit name is part of the CellML spec and available as built-in Sympy quantity
        if unit_name in QuantityStore.CELLML_UNITS and hasattr(units, unit_name):
            self.store[unit_name] = getattr(units, unit_name)
            return self.store[unit_name]

        raise RuntimeError('Cannot find the unit with name "%s"' % unit_name)

    def summarise_units(self, expr: sympy.Expr):
        """Given a Sympy expression, will get all the Quantity objects in the expression and
        collect them together to give a single Sympy expression of the units
        """
        # If this expression is a Quantity itself
        if isinstance(expr, units.Quantity):
            # Return it as it is
            return expr

        # Don't descend into Derivative expression (there are no units within!)
        if isinstance(expr, sympy.Derivative):
            return None

        # If this is an exponential statement and we can evaluate it to a number
        # (i.e. all the units inside the exponential have dropped-out)
        if isinstance(expr, sympy.exp) and isinstance(expr.evalf(), sympy.Number):
            # We say this expression is dimensionless
            return self.get_quantity('dimensionless')

        # Units are always part of a Multiplicative expression
        # TODO: check if all other units are always multiplicative!
        if isinstance(expr, sympy.Mul):
            # We only keep quantities
            keep: List[units.Quantity] = []
            # For each of the multiplication arguments that contain Quantity atoms
            for arg in [x for x in expr.args if x.atoms(units.Quantity)]:
                # Descend into the argument and get the Quantities within
                should_keep = self.summarise_units(arg)
                # Keep anything we find (check - we always should!)
                if should_keep:
                    keep.append(should_keep)
            # Perform the same mathematical operation, but on units only (no symbols or numbers)
            return expr.func(*keep)

        # Otherwise, descend into the expression tree
        return expr.func(*[self.summarise_units(x) for x in expr.args])

    @staticmethod
    def is_equal(quantity_1, quantity_2):
        """Converts one quantity into another to see if they are equal

        :param quantity_1: a Sympy Quantity instance
        :param quantity_2: a Sympy Quantity instance
        """
        return units.convert_to(quantity_1, quantity_2) == quantity_2

    @staticmethod
    def get_conversion_factor(from_unit, to_unit):
        return units.convert_to(from_unit, to_unit).n()

    def _sympify(self, string):
        """Takes a string containing a Sympy expression and evaluates it using the required context
        for handling our units
        """
        # logging.info('sympy.sympify(%s)', string)
        return sympy.sympify(string, locals=self.sympify_context)

    def _make_cellml_quantity(self, name):
        """Will use the CellML unit definition and construct a new Quantity object for that unit
        """
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
                expr = '(%s * %s)' % (unit_element['prefix'], expr)

            if 'exponent' in unit_element:
                expr = '((%s)**%s)' % (expr, unit_element['exponent'])
                # Exponent of the <unit> also affects the dimension
                dimension = '(%s)**%s' % (dimension, unit_element['exponent'])

            # Collect/add this particular <unit> definition
            full_unit_expr.append(expr)
            full_dimension.append('(%s)' % dimension)

        # Have sympy evaluate the collected unit dimension and expression
        full_dimension = self._sympify('*'.join(full_dimension))
        full_unit_expr = self._sympify('*'.join(full_unit_expr))

        # Create and return the Quantity object for this CellML <units>
        quantity = units.Quantity(name, full_dimension, full_unit_expr)
        logging.debug('%s=Quantity("%s", %s, %s)',
                      name, name, full_dimension.args[0], full_unit_expr)
        return quantity
