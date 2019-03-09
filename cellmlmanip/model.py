"""
Classes representing a CellML model and its components
"""
import logging
from collections import OrderedDict
from io import StringIO
from typing import Dict, List, Union

import networkx as nx
import rdflib
import sympy

from cellmlmanip.rdf import create_rdf_node
from cellmlmanip.units import UnitStore


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


# Delimiter for variables name in Sympy expressions: <component><delimiter><name>
SYMPY_SYMBOL_DELIMITER = '$'


class DummyData(object):
    """Holds information about a Dummy placeholder in set of Sympy equations. Dummy symbols are
    used to represent variables from the CellMl model or as placeholders for numbers."""

    def __init__(self, name, units, dummy, initial_value=None,
                 public_interface=None, private_interface=None, number=None,
                 **kwargs):
        # Attributes from the <variable> tag in CellML
        self.name = name
        self.units = units
        self.initial_value = initial_value
        self.public_interface = public_interface
        self.private_interface = private_interface
        self.cmeta_id = kwargs.get('cmeta_id', None)

        # The instance of sympy.Dummy representing this variable in equations
        assert isinstance(dummy, sympy.Dummy)
        self.dummy = dummy

        # The sympy.Dummy assigned in place of this variable (via a connection)
        # If they don't have a public or private 'in' interface
        if private_interface != 'in' and public_interface != 'in':
            # This variable will not connect to anything
            self.assignment = self.dummy
        else:
            # This variable will be connected to another variable
            self.assignment = None

        self.type = None

        if number is not None:
            assert isinstance(number, sympy.Number)
            self.type = 'number'

        self.number = number

    def __str__(self) -> str:
        return '%s(%s)' % (
            type(self).__name__,
            ', '.join('%s=%s' % item for item in vars(self).items() if item[1] is not None)
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
        self.rdf: rdflib.Graph = rdflib.Graph()
        self.graph: nx.DiGraph = None

        self.dummy_data: Dict[sympy.Dummy, DummyData] = OrderedDict()
        self.name_to_symbol: Dict[str, sympy.Dummy] = dict()

        self.equations: List[sympy.Eq] = list()

    def add_unit(self, units_name: str, unit_attributes: List[Dict] = None, base_units=False):
        """Adds information about <units> in <model>
        """
        assert not (unit_attributes and base_units), 'Cannot define base unit with unit attributes'
        if base_units:
            self.units.add_base_unit(units_name)
        else:
            self.units.add_custom_unit(units_name, unit_attributes)

    def add_equation(self, equation: sympy.Eq):
        """Add an equation to this model. Equation must be a Sympy equality"""
        assert isinstance(equation, sympy.Eq), 'Equation expression must be equality'
        self.equations.append(equation)

    def add_number(self, dummy: sympy.Dummy, attributes: Dict):
        """Add metadata about a dummy symbol that represents a number in equations"""
        assert isinstance(dummy, sympy.Dummy)
        assert dummy not in self.dummy_data
        assert 'sympy.Number' in attributes
        assert isinstance(attributes['sympy.Number'], sympy.Number)

        name = '%sNum%s%s' % (SYMPY_SYMBOL_DELIMITER, dummy.dummy_index, SYMPY_SYMBOL_DELIMITER)

        self.dummy_data[dummy] = DummyData(name,
                                           self.units.get_quantity(attributes['cellml:units']),
                                           dummy,
                                           number=attributes['sympy.Number'])

    def add_variable(self, *, name, units, initial_value=None,
                     public_interface=None, private_interface=None, **kwargs):
        """Add information about a variable that represents a symbol in equations. Returns the
        sympy.Dummy created by the model to represent the variable in equations"""
        assert name not in self.name_to_symbol, 'Variable %s already exists' % name

        variable = DummyData(name=name,
                             units=self.units.get_quantity(units),
                             dummy=sympy.Dummy(name),
                             initial_value=initial_value,
                             public_interface=public_interface,
                             private_interface=private_interface,
                             **kwargs)
        self.dummy_data[variable.dummy] = variable
        self.name_to_symbol[variable.name] = variable.dummy
        return variable.dummy

    def get_dummy_data(self, name_or_instance: Union[str, sympy.Dummy]) -> DummyData:
        """Look up dummy data for given symbol. Accepts a string name or instance of sympy.Dummy"""
        if isinstance(name_or_instance, str):
            return self.dummy_data[self.name_to_symbol[name_or_instance]]

        assert isinstance(name_or_instance, sympy.Dummy)
        return self.dummy_data[name_or_instance]

    def check_dummy_instances(self):
        """Check that every symbol in list of equations has a metadata entry. Returns two lists:
        ({set of dummy instances}, {dummy instances without metadata entry}"""

        # collect list of all dummy instances
        dummy_instances = set()
        not_found = set()

        for index, equation in enumerate(self.equations):
            atoms = equation.atoms()
            for atom in atoms:
                # we allow NegativeOne because it is used by sympy to perform division
                # and represent negative numbers
                if not isinstance(atom, sympy.numbers.NegativeOne):
                    dummy_instances.add(atom)
                    if atom not in self.dummy_data:
                        logger.critical('%s not found for eq. %s (%s)' % (atom, index, equation))
                        not_found.add(atom)

        return dummy_instances, not_found

    def connect_variables(self, source_variable: str, target_variable: str):
        """Given the source and target component and variable, create a connection by assigning
        the symbol from the source to the target. If units are not the same, it will add an equation
        to the target component reflecting the relationship. If a symbol has not been assigned to
        the source variable, then return False.

        :param source_variable: source variable name
        :param target_variable: target variable name
        """
        logger.debug('Model.connect_variables(%s ⟶ %s)', source_variable, target_variable)

        source_variable = self.get_dummy_data(source_variable)
        target_variable = self.get_dummy_data(target_variable)

        # If the source variable has already been assigned a final symbol
        if source_variable.assignment:

            if target_variable.assignment:
                raise ValueError('Target already assigned to %s before assignment to %s' %
                                 (target_variable.assignment, source_variable.assignment))

            # If source/target variable is in the same unit
            if source_variable.units == target_variable.units:
                # Direct substitution is possible
                target_variable.assignment = source_variable.assignment
                # everywhere the target variable is used, replace with source variable
                for index, equation in enumerate(self.equations):
                    self.equations[index] = equation.xreplace(
                        {target_variable.dummy: source_variable.assignment}
                    )
            else:
                # Requires a conversion, so we add an equation assigning the target dummy variable
                # to the source variable
                self.equations.append(
                    # TODO: do unit conversion here?
                    sympy.Eq(target_variable.dummy, source_variable.assignment)
                )
                logger.info('Connection req. unit conversion: %s', self.equations[-1])

                # The assigned symbol for this variable is itself
                target_variable.assignment = target_variable.dummy
            logger.debug('Updated target: %s', target_variable)
            return True
        # The source variable has not been assigned a symbol, so we can't make this connection
        return False

    def add_rdf(self, rdf: str):
        """ Takes RDF string and stores it in an RDFlib.Graph for the model. Can be called
        repeatedly.
        """
        self.rdf.parse(StringIO(rdf))

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
                state_variable = self.get_dummy_data(state_symbol)
                Model._set_variable_type(state_variable, 'state')

                # Get the free symbol and update the variable information
                free_symbol = set(lhs_symbol.canonical_variables.keys()).pop()
                free_variable = self.get_dummy_data(free_symbol)
                Model._set_variable_type(free_variable, 'free')

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
                    if variable.type in ['state', 'free']:
                        graph.add_node(rhs_symbol, equation=None, variable_type=variable.type)
                        graph.add_edge(rhs_symbol, lhs_symbol)
                    else:
                        # TODO: <variable> with initial_value is a parameter & variable without
                        # any variables on RHS is a parameter
                        # The variable is a constant or parameter of the model
                        # TODO: Can we tell the difference between a parameter and a constant?
                        # TODO: change to "self." once collect units is in Model class
                        rhs_variable = self.get_dummy_data(rhs_symbol)
                        if rhs_variable.number is None:
                            variable.type = 'parameter'
                            unit = rhs_variable.units
                            dummy = sympy.Dummy(str(variable.initial_value))
                            number = sympy.Float(variable.initial_value)
                            self.add_number(dummy, {'cellml:units': str(unit),
                                                    'sympy.Number': number})
                            graph.add_node(rhs_symbol,
                                           equation=sympy.Eq(rhs_symbol, dummy),
                                           variable_type='parameter')
                            graph.add_edge(rhs_symbol, lhs_symbol)

        for node in graph.nodes:
            if not node.is_Derivative:
                variable = self.find_variable({'dummy': node})
                assert len(variable) == 1
                variable = variable.pop()
                for key in ['cmeta_id', 'name', 'units']:
                    if getattr(variable, key):
                        graph.nodes[node][key] = getattr(variable, key)
                if graph.nodes[node].get('variable_type', '') == 'state':
                    if variable.initial_value:
                        graph.nodes[node]['initial_value'] = sympy.Float(variable.initial_value)

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
                    dummy_data = self.get_dummy_data(dummy)
                    if dummy_data.number is not None:
                        subs_dict[dummy] = dummy_data.number

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
        elif not symbols:
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
        if not variable.type:
            variable.type = variable_type
        else:
            logger.warning('Variable %s already has type=="%s". Skip setting "%s"',
                           variable.dummy, variable.type, variable_type)

    def get_symbols(self, expr):
        """Returns the symbols in an expression"""
        symbols = set()
        if expr.is_Derivative or (expr.is_Dummy and not self.get_dummy_data(expr).number):
            symbols.add(expr)
        else:
            for arg in expr.args:
                symbols |= self.get_symbols(arg)
        return symbols

    def find_variable(self, search_dict):
        """Searches for variables by attributes in their DummyData record. Pass a dictionary of
        {name: value}, where name is one the attributes in DummyData """
        matches = []

        # for each component in this model
        for variable in self.dummy_data.values():
            matched = True
            for search_key, search_value in search_dict.items():
                if not (hasattr(variable, search_key) and
                        getattr(variable, search_key) == search_value):
                    matched = False
                    break
            if matched:
                matches.append(variable)
        return matches
