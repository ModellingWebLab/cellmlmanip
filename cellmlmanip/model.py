"""
Classes representing a CellML model and its components
"""
import logging
from collections import OrderedDict, defaultdict, deque
from io import StringIO
from typing import Dict, List, Set, Tuple

import networkx as nx
import rdflib
import sympy

from cellmlmanip.rdf import create_rdf_node
from cellmlmanip.units import UnitStore


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# Delimiter for variables name in Sympy expressions: <component><delimiter><name>
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

        self.variables: Dict[str, Dict] = OrderedDict()
        self.equations: List[sympy.Eq] = []
        self.numbers: Dict[sympy.Dummy, Dict] = {}

    def __str__(self):
        """Pretty-print object"""
        return 'Component(\n  name: %s\n  equations: %s\n  variables: %s\n  numbers: %s\n)' % (
            self.name, self.equations, self.variables, self.numbers
        )

    def set_parent(self, parent_name):
        """Sets the parent of this component"""
        if self.parent:
            raise ValueError('Parent of component %s already %s. Cannot set %s!' % (self.name,
                                                                                    self.parent,
                                                                                    parent_name))
        self.parent = parent_name

    def add_sibling(self, sibling_name):
        """Adds a sibling for this component"""
        if sibling_name in self.siblings:
            raise ValueError('Sibling component %s already added!' % sibling_name)
        self.siblings.add(sibling_name)

    def add_encapsulated(self, encapsulated_name):
        """Adds an encapsulated component to this component"""
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
                # only allow first-order derivatives
                assert len(derivative.variables) == 1
                dummy_symbols.add(derivative.variables[0])

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
                self.variables[variable_name]['dummy'] = symbol
            else:
                raise KeyError('Variable "%s" in component "%s" could not be found.' %
                               (variable_name, self.name))

    def collect_units(self, expr):
        """Descends into the given Sympy expression and returns the appropriate units (if any)"""
        logger.debug('collect_units(%s)', expr)

        expr_unit_info = self._find_unit_info(expr)
        if expr_unit_info:
            expr_unit_info = self.model.units.get_quantity(expr_unit_info)
            unit_info = {expr: expr_unit_info}
            logger.debug('Found unit for "%s" ⟶ "%s"', expr.name, expr_unit_info)
        else:
            unit_info: Dict = {}

        # Traverse down the arguments of this expression, and collect their unit information
        for args in expr.args:
            unit_info.update(self.collect_units(args))

        return unit_info

    def add_units_to_all_equations(self):
        """Iterates over equations in this component and returns the unit and number information
        for each sympy.Dummy symbol"""
        # TODO: rename method to get_dummy_info()
        component_dummy_info = dict()

        # for each equation in this component
        for equation in self.equations:
            # get the unit information for each dummy symbol
            unit_info = self.collect_units(equation)

            # bundle the information in a dict of dicts
            unit_info = {k: {'unit': v} for k, v in unit_info.items()}

            # for each symbol, add number information if necessary
            for dummy in unit_info.keys():
                if dummy in self.numbers:
                    unit_info[dummy]['number'] = self.numbers[dummy]['sympy.Number']

            # add each dummy symbol in this equation to the collected dummy info for the component
            for dummy, dummy_info in unit_info.items():
                if dummy in component_dummy_info:
                    # TODO: check for conflicting information
                    pass
                else:
                    component_dummy_info[dummy] = dummy_info

        # information about all dummy symbols used in equations in this component
        return component_dummy_info

    def _find_unit_info(self, expr):
        """Takes an expression (part of an equation) and searches in the component variables
        and numbers to get the associated information about this expression/symbol

        :param expr: a Sympy expression. The interesting ones are leaf nodes of the expression
        """
        variable_with_assigned = self.find_variable({'assignment': expr})
        if variable_with_assigned:
            return variable_with_assigned[0]['units']

        variable_with_dummy = self.find_variable({'dummy': expr})
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
            raise ValueError('Search criteria cannot be empty.')

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


class Variable(object):
    """Represents a <variable> in a CellML model. A component may contain any number of <variable>
    elements, which define variables that may be mathematically related in the equation blocks
    contained in the component. """
    def __init__(self, *, name, units, dummy: sympy.Symbol=None, initial_value=None,
                 public_interface=None, private_interface=None, **kwargs):
        self.name = name
        self.units = units
        self.initial_value = initial_value
        self.public_interface = public_interface
        self.private_interface = private_interface
        self.dummy = dummy
        self.cmeta_id = kwargs.get('cmeta:id', None)
        self._component_name = kwargs.get('_component_name', None)
        self._original_name = kwargs.get('_original_name', None)

    def __str__(self) -> str:
        return '%s(%s)' % (
            type(self).__name__,
            ', '.join('%s=%s' % item for item in vars(self).items() if item[1])
        )


class NumberWrapper(object):
    def __init__(self, dummy, units, number):
        self.dummy = dummy
        self.units = units
        self.number = number

    def __str__(self) -> str:
        return '%s(%s)' % (
            type(self).__name__,
            ', '.join('%s=%s' % item for item in vars(self).items() if item[1])
        )


class Model(object):
    """Holds all information about a CellML model and exposes it for manipulation (intention!)
    """
    def __init__(self, name: str) -> None:
        """Create a new instance of Model object
        :param name: the name of the model e.g. from <model name="">
        """
        self.name: str = name
        self.units: 'UnitStore' = UnitStore(model=self)
        self.components: Dict[str, Component] = OrderedDict()
        self.connections: List[Tuple] = []
        self.rdf: rdflib.Graph = rdflib.Graph()
        self.dummy_info: Dict[Dict] = defaultdict(dict)
        self.graph: nx.DiGraph = None
        self.variables_x: Dict[str, Variable] = OrderedDict()
        self.equations_x: List[sympy.Eq] = list()
        self.numbers_x: Dict[sympy.Dummy, NumberWrapper] = dict()

    def __str__(self):
        """Pretty-print each of the components in this model"""
        return '\n'.join([str(v) for v in self.components.values()])

    def add_unit(self, units_name: str, unit_attributes: List[Dict] = None, base_units=False):
        """Adds information about <units> in <model>
        """
        assert not (unit_attributes and base_units), 'Cannot define base unit with unit attributes'
        if base_units:
            self.units.add_base_unit(units_name)
        else:
            self.units.add_custom_unit(units_name, unit_attributes)

    def add_equation(self, equation: sympy.Eq):
        assert isinstance(equation, sympy.Eq)
        self.equations_x.append(equation)

    def add_number(self, dummy: sympy.Dummy, attributes: Dict):
        assert 'sympy.Number' in attributes
        assert isinstance(dummy, sympy.Dummy)
        assert dummy not in self.numbers_x
        self.numbers_x[dummy] = NumberWrapper(dummy,
                                              attributes['cellml:units'],
                                              attributes['sympy.Number'])

    def add_component(self, component: Component):
        """Adds name to list of <component>s in the <model>
        """
        if component.name in self.components:
            raise ValueError('Component "%s" already exists. Check CellML.' % component.name)

        self.components[component.name] = component

    def add_variable(self, variable: Variable):
        assert variable.name not in self.variables, 'Variable %s already exists' % variable.name
        assert variable.dummy is None, 'Variable must not have a Sympy Symbol registered'

        variable.dummy = sympy.Dummy(variable.name)
        self.variables_x[variable.name] = variable
        return variable.dummy

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
    def _is_not_assigned(variable_attributes: Dict):
        """Checks whether a variable gets its value from another component. There are two
        possibilities. Either it has:
            (i) public_interface="out" or
            (ii) private_interface="out" without public_interface="in"
        """
        return (
                ('private_interface', 'in') not in variable_attributes.items() and
                ('public_interface', 'in') not in variable_attributes.items()
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
            if 'dummy' not in variable:
                variable['dummy'] = sympy.Dummy(
                    variable['name']
                )

            # If this variable does not get its value from another component
            if Model._is_not_assigned(variable):
                # The assigned variable to this variable *is* the dummy symbol (no assignment)
                variable['assignment'] = variable['dummy']

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
            logger.debug('Try to connect %s and %s', *connection)
            success = self._connect(connection)
            if not success:
                logger.debug('Cannot connect %s (source does not have assignment).', connection)
                connections_to_process.append(connection)
            # TODO: track looping and break if we can't exit

        # All connections have been made, now substitute the original dummy symbol with new dummy
        for component, variable in self.variables:
            for index, equation in enumerate(component.equations):
                # If this variable is used in an equation, it will have been set a sympy.Dummy
                if 'dummy' in variable:
                    # If the variable has been assigned [a new dummy symbol]
                    if 'assignment' in variable:
                        # Replace the original dummy with the assign dummy symbol
                        component.equations[index] = equation.xreplace(
                            {variable['dummy']: variable['assignment']}
                        )
                    else:
                        variable['assignment'] = variable['dummy']
                # The variable does not have a sympy.Dummy variable set - why??
                else:
                    logger.warning('Variable (%s) in component (%s) not assigned a dummy',
                                   variable['name'],
                                   component.name)

    def add_units_to_equations(self):
        """Inserts unit attributes associated with symbols into equations.

        WARNING: This method replaces/mutates the existing equations i.e. the original equations
        are lost.
        """
        # TODO: rename this method to "collect_dummy_information"

        # this method should only be called once
        if self.dummy_info:
            raise ValueError('Model has already has dummy information.')

        # For each equation in this component
        for component in self.components.values():
            component_dummy_info = component.add_units_to_all_equations()
            for key, value in component_dummy_info.items():
                if key in self.dummy_info:
                    # TODO: check it's the same information
                    pass
                else:
                    self.dummy_info[key] = value

    def check_left_right_units_equal(self, equality: sympy.Eq):
        """Given a Sympy Equality expression, checks that the LHS and RHS have the same units"""
        rhs: sympy.Expr = equality.rhs
        lhs: sympy.Expr = equality.lhs

        lhs_units = self.units.summarise_units(lhs)
        rhs_units = self.units.summarise_units(rhs)

        assert self.units.is_unit_equal(rhs_units, lhs_units), 'Units %s %s != %s %s' % (
            lhs_units, self.units.ureg.get_base_units(lhs_units),
            rhs_units, self.units.ureg.get_base_units(rhs_units)
        )

    def get_equation_graph(self, refresh=False) -> nx.DiGraph:
        """Returns an ordered list of equations for the model"""
        # TODO: Set the parameters of the model (parameters rather than use initial values)

        # if we already have generated the equation graph
        if self.graph and not refresh:
            # return the cached object
            return self.graph

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
                state_variable = self.find_variable({'dummy': state_symbol})
                assert len(state_variable) == 1
                Model._set_variable_type(state_variable[0], 'state')

                # Get the free symbol and update the variable information
                free_symbol = set(lhs_symbol.canonical_variables.keys()).pop()
                free_variable = self.find_variable({'dummy': free_symbol})
                assert len(free_variable) == 1
                Model._set_variable_type(free_variable[0], 'free')

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
                    variable = self.find_variable({'dummy': rhs_symbol})
                    assert len(variable) == 1
                    variable = variable[0]

                    # If the variable is a state or free variable of a derivative
                    if 'type' in variable and variable['type'] in ['state', 'free']:
                        graph.add_node(rhs_symbol, equation=None, variable_type=variable['type'])
                        graph.add_edge(rhs_symbol, lhs_symbol)
                    else:
                        # TODO: <variable> with initial_value is a parameter & variable without
                        # any variables on RHS is a parameter
                        # The variable is a constant or parameter of the model
                        variable['type'] = 'parameter'
                        # TODO: Can we tell the difference between a parameter and a constant?
                        # TODO: change to "self." once collect units is in Model class
                        unit = next(iter(self.components.values())).collect_units(rhs_symbol)
                        unit = unit[rhs_symbol]
                        number = sympy.Dummy(str(variable['initial_value']))
                        self.dummy_info[number] = {
                            'unit': unit,
                            'number': sympy.Float(variable['initial_value'])
                        }
                        graph.add_node(rhs_symbol,
                                       equation=sympy.Eq(rhs_symbol, number),
                                       variable_type='parameter')
                        graph.add_edge(rhs_symbol, lhs_symbol)

        for node in graph.nodes:
            if not node.is_Derivative:
                variable = self.find_variable({'dummy': node})
                assert len(variable) == 1
                variable = variable.pop()
                for key in ['cmeta_id', 'name', 'units']:
                    if key in variable:
                        graph.nodes[node][key] = variable[key]
                if graph.nodes[node].get('variable_type', '') == 'state':
                    if 'initial_value' in variable:
                        graph.nodes[node]['initial_value'] = sympy.Float(variable['initial_value'])

        # TODO: the replacement of dummy-number placeholders can happen above?
        # for each node in the graph
        for node in graph.nodes:
            # if an equation exists for this node
            equation = graph.nodes[node]['equation']
            if equation is not None:
                # get all the dummy symbols on the RHS
                dummies = equation.rhs.atoms(sympy.Dummy)

                # get any dummy-numbers
                subs_dict = {}
                for dummy in dummies:
                    if 'number' in self.dummy_info[dummy]:
                        subs_dict[dummy] = self.dummy_info[dummy]['number']

                # if there are any dummy-numbers on the rhs
                if subs_dict:
                    # replace the equation with a new equation with rhs subbed with real numbers
                    graph.nodes[node]['equation'] = sympy.Eq(equation.lhs,
                                                             equation.rhs.subs(subs_dict))

        self.graph = graph
        return graph

    def get_equations_for(self, symbols):
        graph = self.get_equation_graph()
        sorted_symbols = nx.lexicographical_topological_sort(graph, key=str)

        # Create set of symbols for which we require equations
        required_symbols = set()
        for output in symbols:
            required_symbols.add(output)
            required_symbols.update(nx.ancestors(graph, output))

        eqs = []
        for symbol in sorted_symbols:
            # Ignore symbols we don't need
            if symbol not in required_symbols:
                continue

            # Get equation
            eq = graph.nodes[symbol]['equation']

            # Skip symbols that are not set with an equation
            if eq is None:
                continue

            eqs.append(eq)

        return eqs

    def get_derivative_symbols(self):
        """
        Returns a list of derivative symbols found in the given model graph.
        """
        return [v for v in self.graph if isinstance(v, sympy.Derivative)]

    def get_state_symbols(self):
        """
        Returns a list of state variables found in the given model graph.
        """
        return [v.args[0] for v in self.get_derivative_symbols()]

    def get_free_variable_symbol(self):
        """
        Returns the free variable of the given model graph.
        """
        for v in self.graph:
            if self.graph.nodes[v].get('variable_type', '') == 'free':
                return v

        # This should be unreachable
        raise ValueError('No free variable set in model.')  # pragma: no cover

    def get_symbol_by_cmeta_id(self, cmeta_id):
        """
        Searches the given graph and returns the symbol for the variable with the
        given cmeta_id.
        """
        # TODO: Either add an argument to allow derivative symbols to be fetched, or
        #      create a separate method for them.
        for v in self.graph:
            if self.graph.nodes[v].get('cmeta_id', '') == cmeta_id:
                return v

        raise KeyError('No variable with cmeta id "%s" found.' % str(cmeta_id))

    def get_symbol_by_ontology_term(self, namespace_uri, local_name):
        """
        Searches the RDF graph for a variable annotated with the given
        ``{namespace_uri}local_name`` and returns its symbol.

        Specifically, this method searches for a unique variable annotated with
        predicate ``http://biomodels.net/biology-qualifiers/is`` and the object
        specified by ``{namespace_uri}local_name``.

        Will raise a ``KeyError`` if no variable with the given annotation is
        found, and a ``ValueError`` if more than one variable with the given
        annotation is found.
        """
        symbols = self._get_symbols_by_rdf(
            ('http://biomodels.net/biology-qualifiers/', 'is'),
            (namespace_uri, local_name))
        if len(symbols) == 1:
            return symbols[0]
        elif len(symbols) == 0:
            raise KeyError(
                'No variable annotated with {' + namespace_uri + '}'
                + local_name + ' found.')
        else:
            raise ValueError(
                'Multiple variables annotated with {' + namespace_uri + '}'
                + local_name)

    def _get_symbols_by_rdf(self, predicate, object_=None):
        """
        Searches the RDF graph for variables annotated with the given predicate
        and object (e.g. "is oxmeta:time") and returns the associated symbols.

        Both ``predicate`` and ``object_`` (if given) must be
        ``(namespace, local_name)`` tuples.
        """
        # Convert property and value to RDF nodes
        # TODO: Eventually a different form for predicate and object may be
        #       accepted.
        assert len(predicate) == 2
        predicate = create_rdf_node(*predicate)
        if object_ is not None:
            assert len(object_) == 2
            object_ = create_rdf_node(*object_)

        # Find symbols
        symbols = []
        for result in self.rdf.subjects(predicate, object_):
            assert isinstance(result, rdflib.URIRef), 'Non-resource annotated.'

            # Get cmeta id from result uri
            uri = str(result)
            if uri[0] != '#':
                # TODO This should eventually be implemented
                raise NotImplementedError(
                    'Non-local annotations are not supported.')
            symbols.append(self.get_symbol_by_cmeta_id(uri[1:]))

        return symbols

    def get_value(self, symbol):
        """
        Returns the evaluated value of the given symbol's RHS.
        """
        # Find RHS
        rhs = self.graph.nodes[symbol]['equation'].rhs

        # Evaluate and return
        return float(rhs.evalf())

    def get_initial_value(self, symbol):
        """
        Returns the initial value of the given symbol
        :param symbol: Sympy Dummy object of required symbol
        :return: float of initial value
        """
        return float(self.graph.nodes[symbol]['initial_value'])

    @staticmethod
    def _set_variable_type(variable, variable_type):
        if 'type' not in variable:
            variable['type'] = variable_type
        else:
            logger.warning('Variable %s already has type=="%s". Skip setting "%s"',
                           variable['dummy'], variable['type'], variable_type)

    def get_symbols(self, expr):
        """Returns the symbols in an expression"""
        symbols = set()
        if expr.is_Derivative or (expr.is_Dummy and 'number' not in self.dummy_info[expr]):
            symbols.add(expr)
        else:
            for arg in expr.args:
                symbols |= self.get_symbols(arg)
        return symbols

    def _get_connection_endpoints(self, connection):
        """Pull out the variable dict of the component for the two endpoints of the connection

        :param connection: single connection tuple as created by the parser
        """
        ((component_1, variable_1), (component_2, variable_2)) = connection
        variable_1_attributes = self.components[component_1].variables[variable_1]
        variable_2_attributes = self.components[component_2].variables[variable_2]
        return component_1, variable_1_attributes, component_2, variable_2_attributes

    def _connect(self, connection):
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
        comp_1, var_1, comp_2, var_2 = self._get_connection_endpoints(connection)

        # keys for lookup
        pub = 'public_interface'
        pri = 'private_interface'

        def _are_siblings(comp_a, comp_b):
            return self.components[comp_a].parent == self.components[comp_b].parent

        def _parent_of(parent_name, child_name):
            return parent_name == self.components[child_name].parent

        def _has_interface(dic, key, val):
            return key in dic and dic[key] == val

        # if the components are siblings (either same parent or top-level)
        if _are_siblings(comp_1, comp_2):
            # they are both connected on their public_interface
            if _has_interface(var_1, pub, 'out') and _has_interface(var_2, pub, 'in'):
                return self._connect_with_direction(comp_1, var_1, comp_2, var_2)
            elif _has_interface(var_1, pub, 'in') and _has_interface(var_2, pub, 'out'):
                return self._connect_with_direction(comp_2, var_2, comp_1, var_1)
        else:
            # determine which component is parent of the other
            if _parent_of(comp_1, comp_2):
                parent_comp, child_comp = comp_1, comp_2
                parent_var, child_var = var_1, var_2
            else:
                parent_comp, child_comp = comp_2, comp_1
                parent_var, child_var = var_2, var_1

            # parent/child components are connected using private/public interface, respectively
            if _has_interface(child_var, pub, 'in') and _has_interface(parent_var, pri, 'out'):
                return self._connect_with_direction(parent_comp, parent_var, child_comp, child_var)
            elif _has_interface(child_var, pub, 'out') and _has_interface(parent_var, pri, 'in'):
                return self._connect_with_direction(child_comp, child_var, parent_comp, parent_var)

        raise ValueError('Cannot determine the source & target for connection %s' % str(connection))

    def _connect_with_direction(self,
                                source_component, source_variable,
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
        logger.debug('    Source: %s ⟶ %s', source_component, source_variable)
        logger.debug('    Target: %s ⟶ %s', target_component, target_variable)

        # If the source variable has already been assigned a final symbol
        if 'assignment' in source_variable:

            if 'assignment' in target_variable:
                raise ValueError('Target already assigned to %s before assignment to %s' %
                                 (target_variable['assignment'], source_variable['assignment']))

            # If source/target variable is in the same unit
            if source_variable['units'] == target_variable['units']:
                # Direct substitution is possible
                target_variable['assignment'] = source_variable['assignment']
            else:
                # Requires a conversion, so we add an equation to the component that assigns the
                # target dummy variable to the source variable (unit conversion handled separately)
                self.components[target_component].equations.append(
                    sympy.Eq(target_variable['dummy'], source_variable['assignment'])
                )
                logger.info('    New target eq: %s ⟶ %s', target_component,
                            self.components[target_component].equations[-1])

                # The assigned symbol for this variable is itself
                target_variable['assignment'] = target_variable['dummy']
            logger.debug('    Updated target: %s ⟶ %s', target_component, target_variable)
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
