"""
Classes representing a CellML model and its components
"""
import logging
from collections import deque
from io import StringIO
from typing import Dict, List, Set, Tuple

import networkx as nx
import rdflib
import sympy

from cellmlmanip.units import UnitStore, UnitDummy

# Delimiter for the name of the Sympy symbol: <component><delimiter><name>
SYMPY_SYMBOL_DELIMITER = '$'


class Component(object):
    """Represents a <component> element in CellML <model>"""

    def __init__(self, name: str, model: 'Model') -> None:
        """Create a new Component object

        :param name: the name of the component, usually from <component name="">
        :param model: a reference to the parent Model
        """
        self.name: str = name
        self.model: 'Model' = model

        # store for <group> relationships
        self.parent: str = None
        self.siblings: Set[str] = set()
        self.encapsulated: Set[str] = set()

        self.variables: Dict[str, Dict] = {}
        self.equations: List[sympy.Eq] = []
        self.numbers: Dict[sympy.Dummy, Dict] = {}

    def __str__(self):
        """Pretty-print object"""
        return "Component(\n  name: %s\n  equations: %s\n  variables: %s\n  numbers: %s\n)" % (
            self.name, self.equations, self.variables, self.numbers
        )

    def set_parent(self, parent_name):
        if self.parent:
            raise ValueError('Parent of component %s already %s. Cannot set %s!' % (self.name,
                                                                                    self.parent,
                                                                                    parent_name))
        self.parent = parent_name

    def add_sibling(self, sibling_name):
        if sibling_name in self.siblings:
            raise ValueError('Sibling component %s already added!' % sibling_name)
        self.siblings.add(sibling_name)

    def add_encapsulated(self, encapsulated_name):
        if encapsulated_name in self.encapsulated:
            raise ValueError('Encapsulated component %s already added!' % encapsulated_name)
        self.encapsulated.add(encapsulated_name)

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
            variable_name = symbol.name.split(SYMPY_SYMBOL_DELIMITER)[1]
            if variable_name in self.variables:
                self.variables[variable_name]['sympy.Dummy'] = symbol
            else:
                raise KeyError('Variable "%s" in component "%s" could not be found.' %
                               (variable_name, self.name))

    def collect_units(self, expr):
        """Descends into the given Sympy expression and returns the appropriate units (if any)"""
        logging.debug("collect_units(%s)", expr)

        expr_unit_info = self._find_unit_info(expr)
        if expr_unit_info:
            expr_unit_info = self.model.units.get_quantity(expr_unit_info)
            unit_info = {expr: expr_unit_info}
            logging.debug('Found unit for "%s" ⟶ "%s"', expr.name, expr_unit_info)
        else:
            unit_info: Dict = {}

        # Traverse down the arguments of this expression, and collect their unit information
        for args in expr.args:
            unit_info.update(self.collect_units(args))

        return unit_info

    def get_equation_with_units(self, equation):
        # extract the unit information from the equation
        unit_info = self.collect_units(equation)

        # Re-purpose the unit_info dict into a substitution dictionary
        # to replace the "symbol-or-number" with "unit * symbol-or-number"
        for key, value in unit_info.items():
            # If this unit information is about a number
            if key in self.numbers:
                # Replace the dummy number with the sympy.Number (generated by the parser)
                key.unit = value
                key.number = self.numbers[key]['sympy.Number']
            # Otherwise, straightforward insertion of unit
            else:
                key.unit = value

    def add_units_to_all_equations(self):
        for index, equation in enumerate(self.equations):
            # Warning! This overwrites the original equation
            # self.equations[index] = self.get_equation_with_units(equation)
            self.get_equation_with_units(equation)

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
        # a variable in a different component (above, we're only looking in current component)
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
            raise ValueError("Search criteria cannot be empty.")

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
        self.units: 'UnitStore' = None
        self.components: Dict[str, Component] = {}
        self.connections: List[Tuple] = []
        self.rdf: rdflib.Graph = rdflib.Graph()

    def __str__(self):
        """Pretty-print each of the components in this model"""
        return '\n'.join([str(v) for v in self.components.values()])

    def add_unit(self, units_elements: dict):
        """Adds information about <units> in <model>
        """
        self.units = UnitStore(units_elements)

    def add_component(self, component: Component):
        """Adds name to list of <component>s in the <model>
        """
        if component.name in self.components:
            raise ValueError("%s already exists. Check CellML." % component.name)

        self.components[component.name] = component

    def add_rdf(self, rdf: str):
        """ Takes RDF string and stores it in an RDFlib.Graph for the model. Can be called
        repeatedly.
        """
        self.rdf.parse(StringIO(rdf))

    @property
    def equations(self):
        """Helper method to iterate over equations in the model (across different components)"""
        for component in self.components.values():
            for equation in component.equations:
                yield equation

    @property
    def variables(self):
        """Iterates over all variables in the model (across different components)"""
        for component in self.components.values():
            for var_attr in component.variables.values():
                yield (component, var_attr)

    @staticmethod
    def __is_not_assigned(variables: Dict):
        """Checks whether a variable gets its value from another component. There are two
        possibilities. Either it has:
            (i) public_interface="out" or
            (ii) private_interface="out" without public_interface="in"
        """
        return (
            ('private_interface', 'in') not in variables.items() and
            ('public_interface', 'in') not in variables.items()
        )

    def make_connections(self):
        """Uses public/private interface attributes of variables and model connections to assign
        target variables to their source
        """

        # At this stage we are going to assign variables to new sympy Dummy placeholders
        # This may change in the future once we're confident that everything is working e.g.
        # manipulate the equations directly

        # First, we assign new sympy.Dummy variables for those CellML <variable>s that are only
        # exposed to other components i.e. do not have an "in" value on public/private_interface

        # For each CellML <variable> in the model
        for component, variable in self.variables:
            # If this variable does not have a sympy.Dummy (e.g. variable declared but not mathml?)
            if 'sympy.Dummy' not in variable:
                variable['sympy.Dummy'] = UnitDummy(
                    component.name + SYMPY_SYMBOL_DELIMITER + variable['name']
                )

            # If this variable does not get its value from another component
            if Model.__is_not_assigned(variable):
                # The assigned variable to this variable *is* the dummy symbol (no assignment)
                variable['assignment'] = variable['sympy.Dummy']

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
            if not success:
                logging.debug('Cannot connect %s (source does not have assignment).', connection)
                connections_to_process.append(connection)
            # TODO: track looping and break if we can't exit

        # All connections have been made, now substitute the original dummy symbol with new dummy
        for component, variable in self.variables:
            for index, equation in enumerate(component.equations):
                # If this variable is used in an equation, it will have been set a sympy.Dummy
                if 'sympy.Dummy' in variable:
                    # If the variable has been assigned [a new dummy symbol]
                    if 'assignment' in variable:
                        # Replace the original dummy with the assign dummy symbol
                        component.equations[index] = equation.xreplace(
                            {variable['sympy.Dummy']: variable['assignment']}
                        )
                    else:
                        variable['assignment'] = variable['sympy.Dummy']
                # The variable does not have a sympy.Dummy variable set - why??
                else:
                    logging.warning('Variable (%s) in component (%s) not assigned a dummy',
                                    variable['name'], component.name)

    def add_units_to_equations(self):
        """Inserts unit attributes associated with symbols into equations.

        WARNING: This method replaces/mutates the existing equations i.e. the original equations
        are lost.
        """
        # For each equation in this component
        for component in self.components.values():
            component.add_units_to_all_equations()

    def check_left_right_units_equal(self, equality: sympy.Eq):
        """Given a Sympy Equality expression, checks that the LHS and RHS have the same units"""
        rhs: sympy.Expr = equality.rhs
        lhs: sympy.Expr = equality.lhs

        if rhs.is_Piecewise:
            for piece, _ in rhs.args:
                self.check_left_right_units_equal(sympy.Eq(lhs, piece))
        else:
            lhs_units = self.units.summarise_units(lhs)
            rhs_units = self.units.summarise_units(rhs)

            assert self.units.is_unit_equal(rhs_units, lhs_units), 'Units %s != %s' % (rhs_units,
                                                                                       lhs_units)

    def get_equation_graph(self) -> nx.DiGraph:
        """Returns an ordered list of equations for the model"""
        # TODO: Set the parameters of the model (parameters rather than use initial values)

        # store symbols, their attributes and their relationships in a directed graph
        graph = nx.DiGraph()

        equation_count = 0

        # for each equation in the model
        for equation in self.equations:
            equation_count += 1

            # Determine LHS. We should only every have one symbol (or derivative)
            lhs_symbol = self.get_symbols(equation.lhs)
            assert len(lhs_symbol) == 1
            lhs_symbol = lhs_symbol.pop()

            # Add the lhs symbol of the equation to the graph
            graph.add_node(lhs_symbol, equation=equation)

            # If LHS is a derivative
            # TODO: should be in perhaps collect_variable_attributes() but after connections made?
            if lhs_symbol.is_Derivative:
                # Get the state symbol and update the variable information
                state_symbol = lhs_symbol.free_symbols.pop()
                state_variable = self.find_variable({'sympy.Dummy': state_symbol})
                assert len(state_variable) == 1
                Model.__set_variable_type(state_variable[0], 'state')

                # Get the free symbol and update the variable information
                free_symbol = set(lhs_symbol.canonical_variables.keys()).pop()
                free_variable = self.find_variable({'sympy.Dummy': free_symbol})
                assert len(free_variable) == 1
                Model.__set_variable_type(free_variable[0], 'free')

        # sanity check none of the lhs have the same hash!
        assert len(graph.nodes) == equation_count

        # sanity check all the lhs are unique in meaning (sympy.Dummy: same name != same hash)
        assert len(set([str(x) for x in graph.nodes])) == equation_count

        # for each equation in the model
        for equation in self.equations:
            # get the lhs symbol
            lhs_symbol = self.get_symbols(equation.lhs).pop()

            # for each of the symbols on the rhs of the equation
            for rhs_symbol in self.get_symbols(equation.rhs):
                # if the symbol maps to a node in the graph
                if rhs_symbol in graph.node:
                    # add the dependency edge
                    graph.add_edge(rhs_symbol, lhs_symbol)
                else:
                    # The symbol does not have a node in the graph, get the variable info
                    variable = self.find_variable({'sympy.Dummy': rhs_symbol})
                    assert len(variable) == 1
                    variable = variable[0]

                    # If the variable is a state or free variable of a derivative
                    if 'type' in variable and variable['type'] in ['state', 'free']:
                        graph.add_node(rhs_symbol, equation=None, variable_type=variable['type'])
                        graph.add_edge(rhs_symbol, lhs_symbol)
                    else:
                        # TODO: <variable> with initial_value is a paramter & variable without any variables on RHS is a parameter
                        # The variable is a constant or parameter of the model
                        variable['type'] = 'parameter'
                        # TODO: Can we tell the difference between a parameter and a constant?
                        # TODO: change to "self." once collect units is in Model class
                        unit = next(iter(self.components.values())).collect_units(rhs_symbol)
                        unit = unit[rhs_symbol]
                        number = UnitDummy(str(variable['initial_value']))
                        number.unit = unit
                        number.number = sympy.Float(variable['initial_value'])
                        graph.add_node(rhs_symbol,
                                       equation=sympy.Eq(rhs_symbol, number),
                                       variable_type='parameter')
                        graph.add_edge(rhs_symbol, lhs_symbol)

        for node in graph.nodes:
            if not node.is_Derivative:
                variable = self.find_variable({'sympy.Dummy': node})
                assert len(variable) == 1
                variable = variable.pop()
                for key in ['cmeta:id', 'name', 'units']:
                    if key in variable:
                        graph.nodes[node][key] = variable[key]
                if graph.nodes[node].get('variable_type', '') == 'state':
                    if 'initial_value' in variable:
                        graph.nodes[node]['initial_value'] = sympy.Float(variable['initial_value'])

        return graph

    @staticmethod
    def __set_variable_type(variable, variable_type):
        if 'type' not in variable:
            variable['type'] = variable_type
        else:
            logging.error("The variable %s has been set a type of '%s'. Skip setting '%s'",
                          variable['sympy.Dummy'], variable['type'], variable_type)

    def get_symbols(self, expr):
        """Returns the symbols in an expression"""
        symbols = set()
        if expr.is_Derivative or (expr.is_Dummy and expr.number is None):
            symbols.add(expr)
        else:
            for arg in expr.args:
                symbols |= self.get_symbols(arg)
        return symbols

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

        Relevant lines from the CellML specification:

            The set of all components immediately encapsulated by the current
            component is the encapsulated set.

            Other components encapsulated by the same parent make up the
            sibling set.

            The interface exposed to the parent component and components in
            the sibling set is defined by the public_interface attribute. The
            private_interface attribute defines the interface exposed to
            components in the encapsulated set. Each interface has three possible
            values: "in", "out", and "none", where "none" indicates the absence
            of an interface.

        :param connection: a single connection tuple, created by the CellML parser
            ((component_1, variable_1), (component_2, variable_2))
        """
        comp_1, var_1, comp_2, var_2 = self.__get_connection_endpoints(connection)

        # keys for lookup
        pub = 'public_interface'
        pri = 'private_interface'

        def __are_siblings(comp_a, comp_b):
            return self.components[comp_a].parent == self.components[comp_b].parent

        def __parent_of(parent_name, child_name):
            return parent_name == self.components[child_name].parent

        def __has_interface(dic, key, val):
            return key in dic and dic[key] == val

        # if the components are siblings (either same parent or top-level)
        if __are_siblings(comp_1, comp_2):
            # they are both connected on their public_interface
            if __has_interface(var_1, pub, 'out') and __has_interface(var_2, pub, 'in'):
                return self.__connect_with_direction(comp_1, var_1, comp_2, var_2)
            elif __has_interface(var_1, pub, 'in') and __has_interface(var_2, pub, 'out'):
                return self.__connect_with_direction(comp_2, var_2, comp_1, var_1)
        else:
            # determine which component is parent of the other
            if __parent_of(comp_1, comp_2):
                parent_comp, child_comp = comp_1, comp_2
                parent_var, child_var = var_1, var_2
            else:
                parent_comp, child_comp = comp_2, comp_1
                parent_var, child_var = var_2, var_1

            # parent/child components are connected using private/public interface, respectively
            if __has_interface(child_var, pub, 'in') and __has_interface(parent_var, pri, 'out'):
                return self.__connect_with_direction(parent_comp, parent_var, child_comp, child_var)
            elif __has_interface(child_var, pub, 'out') and __has_interface(parent_var, pri, 'in'):
                return self.__connect_with_direction(child_comp, child_var, parent_comp, parent_var)

        raise ValueError("Cannot determine the source & target for connection %s" % str(connection))

    def __connect_with_direction(self, source_component, source_variable,
                                 target_component, target_variable):
        """Given the source and target component and variable, create a connection by assigning
        the symbol from the source to the target. If units are not the same, it will add an equation
        to the target component reflecting the relationship. If a symbol has not been assigned to
        the source variable, then return False.

        :param source_component: name of source component
        :param source_variable: name of source variable
        :param target_component: name of target component
        :param target_variable: name of target variable
        """
        logging.debug('    Source: %s ⟶ %s', source_component, source_variable)
        logging.debug('    Target: %s ⟶ %s', target_component, target_variable)

        # If the source variable has already been assigned a final symbol
        if 'assignment' in source_variable:

            if 'assignment' in target_variable:
                raise ValueError('Target already assigned to %s before assignment to %s',
                                 target_variable['assignment'],
                                 source_variable['assignment'])

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


