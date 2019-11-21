"""Classes to represent a flattened CellML model and metadata about its variables"""
import logging
from collections import OrderedDict
from io import StringIO

import networkx as nx
import rdflib
import sympy

from cellmlmanip.rdf import create_rdf_node
from cellmlmanip.units import UnitStore


logger = logging.getLogger(__name__)


# Delimiter for variables name in Sympy expressions: <component><delimiter><name>
SYMPY_SYMBOL_DELIMITER = '$'


class MetaDummy(object):
    """
    Holds information about a Dummy placeholder in set of Sympy equations.

    Dummy symbols are used to represent variables from the CellMl model or as placeholders for numbers.
    """
    # TODO Add param list in docstring above

    def __init__(self, name, units, dummy, initial_value=None,
                 public_interface=None, private_interface=None, number=None,
                 order_added=None,
                 **kwargs):
        # Attributes from the <variable> tag in CellML
        self.name = name
        self.units = units
        self.initial_value = initial_value
        self.public_interface = public_interface
        self.private_interface = private_interface
        self.order_added = order_added
        self.cmeta_id = kwargs.get('cmeta_id', None)

        # The instance of sympy.Dummy representing this variable in equations
        assert isinstance(dummy, sympy.Dummy)
        self.dummy = dummy

        # The sympy.Dummy assigned in place of this variable (via a connection)
        # If they don't have a public or private 'in' interface
        if private_interface != 'in' and public_interface != 'in':
            # This variable is assigned to itself
            self.assigned_to = self.dummy
        else:
            # This variable will be connected to another variable
            self.assigned_to = None

        self.type = None

        if number is not None:
            assert isinstance(number, sympy.Number)
            self.type = 'number'

        self.number = number

    @property
    def is_number(self):
        """Indicates whether this dummy instance is used as a placeholder for a number"""
        return self.number is not None

    def __str__(self):
        return '%s(%s)' % (
            type(self).__name__,
            ', '.join('%s=%s' % item for item in vars(self).items() if item[1] is not None)
        )


class Model(object):
    """
    A componentless representation of a CellML model, containing a list of equations, units, and RDF metadata about
    symbols used in those equations.

    The main parts of a Model are 1. a list of sympy equation objects; 2. a collection of named units; and 3. an RDF
    graph that stores further meta data about the model.

    Equations are stored as ``Sympy.Eq`` objects, but with the caveat that all variables and numbers must be specified
    using the ``Sympy.Dummy`` objects returned by :meth:`add_variable()` and :meth:`add_number()`.

    Cellmlmanip does not support algebraic models: the left-hand side every equation in the model must be a variable or
    a derivative.

    :param name: the name of the model e.g. from ``<model name="">``.
    :param cmeta_id: An optional cmeta id, e.g. from ``<model cmeta:id="">``.
    """
    def __init__(self, name, cmeta_id=None):

        self.name = name
        self.cmeta_id = cmeta_id
        self.rdf_identity = rdflib.URIRef('#' + cmeta_id) if cmeta_id else None

        # A list of sympy.Eq equation objects
        self.equations = []

        # A pint UnitStore, mapping unit names to unit objects
        self.units = UnitStore(model=self)

        # An RDF graph containing further meta data
        self.rdf = rdflib.Graph()

        # Maps sympy.Dummy objects to MetaDummy objects
        self.dummy_metadata = OrderedDict()

        # Maps string names to sympy.Dummy objects
        self.name_to_symbol = dict()

        # A cached nx.DiGraph of this model's equations, created by get_equation_graph()
        self.graph = None

    def add_unit(self, name, attributes=None, base_units=False):
        """
        Adds a unit of measurement to this model, with a given ``name`` and list of ``attributes``.

        :param name: A string name.
        :param attributes: An optional list of dictionaries containing unit attributes. See
            :meth:`UnitStore.add_custom_unit()`.
        :base_units: Set to ``True`` to define a new base unit.
        """
        if base_units:
            if attributes:
                raise ValueError('Base units can not be defined with unit attributes.')
            self.units.add_base_unit(name)
        else:
            self.units.add_custom_unit(name, attributes)

    def add_equation(self, equation):
        """
        Adds an equation to this model.

        The left-hand side (LHS) of the equation must be either a variable symbol or a derivative.

        All numbers and variable symbols used in the equation must have been obtained from this model, e.g. via
        :meth:`add_number()`, :meth:`add_variable()`, or :meth:`get_symbol_by_cmeta_id()`.

        :param equation: A ``sympy.Eq`` object.
        """
        assert isinstance(equation, sympy.Eq), 'The argument `equation` must be a sympy.Eq.'
        self.equations.append(equation)

        # Invalidate the cached graph
        self.graph = None

    def add_number(self, *, number, units, dummy=None):
        """
        Add metadata about a dummy symbol that represents a number in equations.

        :param number: A ``sympy.Number``.
        :param units: A string unit representation.
        :param dummy: An optional ``sympy.Dummy``.

        :return: The ``sympy.Dummy`` object used to represent this number.
        """
        assert isinstance(number, sympy.Number), 'The argument `number` must be a sympy.Number.'

        # Create a dummy object if necessary
        if not dummy:
            dummy = sympy.Dummy(str(number))
        else:
            assert dummy not in self.dummy_metadata

        name = '%sNum%s%s' % (SYMPY_SYMBOL_DELIMITER, dummy.dummy_index, SYMPY_SYMBOL_DELIMITER)

        # store metadata information about the number
        self.dummy_metadata[dummy] = MetaDummy(name=name,
                                               units=self.units.get_quantity(units),
                                               dummy=dummy,
                                               number=number)

        return self.dummy_metadata[dummy].dummy

    # TODO: Do we need the * here?
    def add_variable(self, *, name, units, initial_value=None, public_interface=None, private_interface=None, **kwargs):
        """
        Add a variable to the model and return a Sympy ``Dummy`` object to represent it.

        :param name: A string name.
        :param units: A string units representation.
        :param initial_value: An optional initial value.
        :param public_interface: An optional public interface specifier (only required when parsing CellML).
        :param private_interface: An optional private interface specifier (only required when parsing CellML).
        :param kwargs: Any further keyword arguments will be passed to the :class:`MetaDummy` constructor.

        :return: The ``sympy.Dummy`` created by the model to represent the variable in equations.
        """
        if name in self.name_to_symbol:
            raise ValueError('Variable %s already exists.' % name)

        dummy = sympy.Dummy(name)

        if initial_value is not None:
            initial_value = float(initial_value)

        variable = MetaDummy(name=name,
                             units=self.units.get_quantity(units),
                             dummy=dummy,
                             initial_value=initial_value,
                             public_interface=public_interface,
                             private_interface=private_interface,
                             order_added=len(self.dummy_metadata),
                             **kwargs)

        self.dummy_metadata[dummy] = variable
        self.name_to_symbol[name] = dummy

        return self.dummy_metadata[dummy].dummy

    def get_meta_dummy(self, name_or_instance):
        """
        Look up dummy data for given symbol.

        :param name_or_instance: A string name or a ``sympy.Dummy``.
        :return: A :class:`MetaDummy`.
        """
        if isinstance(name_or_instance, str):
            return self.dummy_metadata[self.name_to_symbol[name_or_instance]]

        assert isinstance(name_or_instance, sympy.Dummy)
        return self.dummy_metadata[name_or_instance]

    def check_dummy_metadata(self):
        """
        Check that every symbol in list of equations has a metadata entry.

        :return: A set of dummy instances and a set of dummy instances without metadata.
        """
        # collect list of all dummy instances
        dummy_instances = set()
        not_found = set()

        for index, equation in enumerate(self.equations):
            atoms = equation.atoms(sympy.Dummy)
            for atom in atoms:
                dummy_instances.add(atom)
                if atom not in self.dummy_metadata:
                    logger.critical('%s not found for eq. %s (%s)' % (atom, index, equation))
                    not_found.add(atom)

        return dummy_instances, not_found

    def check_cmeta_id(self):
        """Checks that every variable with a cmeta_id is a source variable"""
        is_okay = True
        for variable in self.dummy_metadata.values():
            if not variable.is_number:
                if variable.cmeta_id:
                    if variable.dummy != variable.assigned_to:
                        is_okay = False
                        logger.critical('%s has cmeta id but is assigned to %s',
                                        variable.dummy, variable.assigned_to)
        return is_okay

    def check_dummy_assignment(self):
        """Every non-number dummy symbol in the model should be assigned to itself or a source
        variable. The source variable must be assigned to itself"""
        is_okay = True
        for variable in self.dummy_metadata.values():
            if not variable.is_number:
                # either the variable is assigned to itself
                if variable.dummy == variable.assigned_to:
                    continue

                # or the variable is assigned to a source variable
                source_dummy = self.dummy_metadata[variable.assigned_to]

                # the source dummy must be assigned to itself
                if source_dummy.dummy == source_dummy.assigned_to:
                    continue

                is_okay = False
                logger.critical('%s is assigned to %s, which is assigned to %s',
                                variable.dummy,
                                variable.assigned_to,
                                source_dummy.assigned_to)
        return is_okay

    def check_variables_in_equations(self):
        """Every variable we have should have been used in an equation."""
        pass

    def connect_variables(self, source_name: str, target_name: str):
        """Given the source and target component and variable, create a connection by assigning
        the symbol from the source to the target. If units are not the same, it will add an equation
        to the target component reflecting the relationship. If a symbol has not been assigned to
        the source variable, then return False.

        :param source_name: source variable name
        :param target_name: target variable name
        """
        logger.debug('connect_variables(%s ⟶ %s)', source_name, target_name)

        source = self.get_meta_dummy(source_name)
        target = self.get_meta_dummy(target_name)

        # If the source variable has already been assigned a final symbol
        if source.assigned_to:

            if target.assigned_to:
                raise ValueError('Target already assigned to %s before assignment to %s' %
                                 (target.assigned_to, source.assigned_to))

            # If source/target variable is in the same unit
            if source.units == target.units:
                # Direct substitution is possible
                target.assigned_to = source.assigned_to
                # everywhere the target variable is used, replace with source variable
                for index, equation in enumerate(self.equations):
                    self.equations[index] = equation.xreplace(
                        {target.dummy: source.assigned_to}
                    )
            # Otherwise, this connection requires a conversion
            else:
                # Get the scaling factor required to convert source units to target units
                factor = self.units.convert_to(1 * source.units, target.units).magnitude

                # Dummy to represent this factor in equations, having units for conversion
                factor_dummy = self.add_number(number=sympy.Float(factor),
                                               units=str(target.units / source.units))

                # Add an equations making the connection with the required conversion
                self.equations.append(sympy.Eq(target.dummy, source.assigned_to * factor_dummy))

                logger.info('Connection req. unit conversion: %s', self.equations[-1])

                # The assigned symbol for this variable is itself
                target.assigned_to = target.dummy

            logger.debug('Updated target: %s', target)

            return True

        # The source variable has not been assigned a symbol, so we can't make this connection
        logger.info('The source variable has not been assigned to a symbol '
                    '(i.e. expecting a connection): %s ⟶ %s',
                    target.name, source.name)
        return False

    def add_rdf(self, rdf: str):
        """ Takes an RDF string and stores it in the model's RDF graph. """
        self.rdf.parse(StringIO(rdf))

    def check_left_right_units_equal(self, equality):
        """
        Checks whether the LHS and RHS in a ``sympy.Eq`` have the same units.
        :param equality: A ``sympy.Eq``.
        """
        lhs_units = self.units.summarise_units(equality.lhs)
        rhs_units = self.units.summarise_units(equality.rhs)

        assert self.units.is_unit_equal(rhs_units, lhs_units), 'Units %s %s != %s %s' % (
            lhs_units, self.units.ureg.get_base_units(lhs_units),
            rhs_units, self.units.ureg.get_base_units(rhs_units)
        )

    def get_equation_graph(self, refresh=False):
        """
        Returns a ``networkx.DiGraph`` containing the model equations.

        :param refresh: If set to ``True`` the graph will be regenerated, even if a cached version is available.
        :return: An ``nx.Digraph`` of equations.
        """
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

            # Determine LHS.
            lhs = equation.lhs
            if not (lhs.is_Derivative or (lhs.is_Dummy and not self.get_meta_dummy(lhs).number)):
                raise RuntimeError('DAEs are not supported. All equations must be of form `x = ...` or `dx/dt = ...')

            # Add the lhs symbol of the equation to the graph
            graph.add_node(lhs, equation=equation)

            # If LHS is a derivative
            if lhs.is_Derivative:
                # Get the state symbol and update the variable information
                state_symbol = lhs.free_symbols.pop()
                state_variable = self.get_meta_dummy(state_symbol)
                state_variable.type = 'state'

                # Get the free symbol and update the variable information
                free_symbol = lhs.variables[0]
                free_variable = self.get_meta_dummy(free_symbol)
                free_variable.type = 'free'
            else:
                variable = self.get_meta_dummy(lhs)
                variable.type = None

        # sanity check none of the lhs have the same hash!
        assert len(graph.nodes) == equation_count

        # sanity check all the lhs are unique in meaning (sympy.Dummy: same name != same hash)
        assert len(set([str(x) for x in graph.nodes])) == equation_count

        # for each equation in the model
        for equation in self.equations:
            # Get the LHS (already checked above)
            lhs = equation.lhs

            # for each of the symbols or derivatives on the rhs of the equation
            for rhs in self.find_symbols_and_derivatives(equation.rhs):
                # if the symbol maps to a node in the graph
                if rhs in graph.nodes:
                    # add the dependency edge
                    graph.add_edge(rhs, lhs)
                else:
                    # The symbol does not have a node in the graph, get the variable info
                    variable = self.find_variable({'dummy': rhs})
                    assert len(variable) == 1
                    variable = variable[0]

                    # If the variable is a state or free variable of a derivative
                    if variable.type in ['state', 'free']:
                        graph.add_node(rhs, equation=None, variable_type=variable.type)
                        graph.add_edge(rhs, lhs)
                    else:
                        # if the variable on the right-hand side is a number
                        rhs_variable = self.get_meta_dummy(rhs)
                        if rhs_variable.number is None:
                            # this variable is a parameter - add to graph and connect to lhs
                            variable.type = 'parameter'
                            unit = rhs_variable.units
                            number = sympy.Float(variable.initial_value)
                            dummy = self.add_number(number=number, units=str(unit))
                            graph.add_node(rhs, equation=sympy.Eq(rhs, dummy), variable_type='parameter')
                            graph.add_edge(rhs, lhs)

        # add metadata about each node directly to the graph
        # TODO: necessary? remove?
        for node in graph.nodes:
            if not node.is_Derivative:
                variable = self.find_variable({'dummy': node})
                assert len(variable) == 1
                variable = variable.pop()
                for key in ['cmeta_id', 'name', 'units']:
                    if getattr(variable, key):
                        graph.nodes[node][key] = getattr(variable, key)
                if graph.nodes[node].get('variable_type', '') == 'state':
                    if variable.initial_value is not None:
                        graph.nodes[node]['initial_value'] = sympy.Float(variable.initial_value)
                if variable.type is not None:
                    graph.nodes[node]['variable_type'] = variable.type

        # for each node in the graph
        for node in graph.nodes:
            # if an equation exists for this node
            equation = graph.nodes[node]['equation']
            if equation is not None:
                # get all the dummy symbols on the RHS
                dummies = equation.rhs.atoms(sympy.Dummy)

                # get any dummy symbols which are placeholders for numbers
                subs_dict = {}
                for dummy in dummies:
                    dummy_data = self.get_meta_dummy(dummy)
                    if dummy_data.number is not None:
                        subs_dict[dummy] = dummy_data.number

                # if there are any dummy-numbers on the rhs
                if subs_dict:
                    # replace the equation with one with the rhs subbed with real numbers
                    graph.nodes[node]['equation'] = sympy.Eq(equation.lhs,
                                                             equation.rhs.subs(subs_dict))

        self.graph = graph
        return graph

    def get_equations_for(self, symbols, lexicographical_sort=True, recurse=True):
        """Get all equations for a given collection of symbols.

        Results are sorted in topographical order.
        :param symbols: the symbols to get the equations for
        :param lexicographical_sort: indicates whether the result is sorted in lexicographical order first
        :param recurse: indicates whether to recurse the equation graph, or to return only the top level equations
        """
        graph = self.get_equation_graph()
        if lexicographical_sort:
            sorted_symbols = nx.lexicographical_topological_sort(graph, key=str)
        else:
            sorted_symbols = nx.topological_sort(graph)

        # Create set of symbols for which we require equations
        required_symbols = set()
        for output in symbols:
            required_symbols.add(output)
            if recurse:
                required_symbols.update(nx.ancestors(graph, output))
            else:
                required_symbols.update(graph.pred[output])

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
        """Returns a list of derivative symbols found in the given model graph.
        The list is ordered by appearance in the cellml document.
        """
        derivative_symbols = [v for v in self.graph if isinstance(v, sympy.Derivative)]
        return sorted(derivative_symbols, key=lambda state_var: self.get_meta_dummy(state_var.args[0]).order_added)

    def get_state_symbols(self):
        """Returns a list of state variables found in the given model graph.
        The list is ordered by appearance in the cellml document.
        """
        state_symbols = [v.args[0] for v in self.get_derivative_symbols()]
        return sorted(state_symbols, key=lambda state_var: self.get_meta_dummy(state_var).order_added)

    def get_free_variable_symbol(self):
        """Returns the free variable of the given model graph.
        """
        for v in self.graph:
            if self.graph.nodes[v].get('variable_type', '') == 'free':
                return v

        # This should be unreachable
        raise ValueError('No free variable set in model.')  # pragma: no cover

    # TODO: Do we still need this method?
    def get_symbol_by_cmeta_id(self, cmeta_id):
        """
        Searches the given graph and returns the symbol for the variable with the given cmeta id.

        To get symbols from e.g. an oxmeta ontology term, use :meth:`get_symbol_by_ontology_term()`.
        """
        for v in self.graph:
            if self.graph.nodes[v].get('cmeta_id', '') == cmeta_id:
                return v

        raise KeyError('No variable with cmeta id "%s" found.' % str(cmeta_id))

    def get_symbol_by_ontology_term(self, namespace_uri, local_name):
        """Searches the RDF graph for a variable annotated with the given
        ``{namespace_uri}local_name`` and returns its symbol.

        Specifically, this method searches for a unique variable annotated with
        predicate ``http://biomodels.net/biology-qualifiers/is`` and the object
        specified by ``{namespace_uri}local_name``.

        Will raise a ``KeyError`` if no variable with the given annotation is
        found, and a ``ValueError`` if more than one variable with the given
        annotation is found.
        """
        symbols = self.get_symbols_by_rdf(
            ('http://biomodels.net/biology-qualifiers/', 'is'),
            (namespace_uri, local_name))
        if len(symbols) == 1:
            return symbols[0]
        elif not symbols:
            raise KeyError('No variable annotated with {%s}%s found.' %
                           (namespace_uri, local_name))
        else:
            raise ValueError('Multiple variables annotated with {%s}%s' %
                             (namespace_uri, local_name))

    def get_rdf_annotations(self, subject=None, predicate=None, object_=None):
        """Searches the RDF graph and returns 'triples matching the given parameters'

        :param subject: the subject of the triples returned
        :param predicate: the predicate of the triples returned
        :param object_: the object of the triples returned

        ``subject`` ``predicate`` and ``object_`` are optional, if None then any triple matches
        if all are none, all triples are returned
        ``subject`` ``predicate`` and ``object_`` can be anything valid as input to create_rdf_node
        typically an (NS, local) pair, a string or None"""
        subject = create_rdf_node(subject)
        predicate = create_rdf_node(predicate)
        object_ = create_rdf_node(object_)
        return self.rdf.triples((subject, predicate, object_))

    def get_rdf_value(self, subject, predicate):
        """Get the value of an RDF object connected to ``subject`` by ``predicate``.

        :param subject: the object of the triple returned
        :param predicate: the object of the triple returned

        Note: expects exactly one triple to match and the result to be a literal. It's string value is  returned."""
        triples = list(self.get_rdf_annotations(subject, predicate))
        assert len(triples) == 1
        assert isinstance(triples[0][2], rdflib.Literal)
        value = str(triples[0][2]).strip()  # Could make this cleverer by considering data type if desired
        return value

    def get_symbols_by_rdf(self, predicate, object_=None):
        """Searches the RDF graph for variables annotated with the given predicate and object (e.g. "is oxmeta:time")
        and returns the associated symbols sorted in document order.

        Both ``predicate`` and ``object_`` (if given) must be ``(namespace, local_name)`` tuples or string literals.
        """
        predicate = create_rdf_node(predicate)
        object_ = create_rdf_node(object_)

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

        return sorted(symbols, key=lambda sym: self.get_meta_dummy(sym).order_added)

    def get_ontology_terms_by_symbol(self, symbol, namespace_uri=None):
        """Searches the RDF graph for the annotation ``{namespace_uri}annotation_name``
        for the given symbol and returns ``annotation_name`` and optionally restricted
        to a specific ``{namespace_uri}annotation_name``

        Specifically, this method searches for a ``annotation_name`` for
        subject symbol
        predicate ``http://biomodels.net/biology-qualifiers/is`` and the object
        specified by ``{namespace_uri}annotation_name``
        for a specific ``{namespace_uri}`` if set, otherwise for any namespace_uri

        Will return a list of term names.
        """
        ontology_terms = []
        cmeta_id = self.graph.nodes[symbol].get('cmeta_id', None)
        if cmeta_id:
            predicate = ('http://biomodels.net/biology-qualifiers/', 'is')
            predicate = create_rdf_node(predicate)
            for delimeter in ('#', '/'):  # Look for terms using either possible namespace delimiter
                subject = rdflib.term.URIRef(delimeter + cmeta_id)
                for object in self.rdf.objects(subject, predicate):
                    # We are only interested in annotation within the namespace
                    if namespace_uri is None or str(object).startswith(namespace_uri):
                        uri_parts = str(object).split(delimeter)
                        ontology_terms.append(uri_parts[-1])
        return ontology_terms

    def has_ontology_annotation(self, symbol, namespace_uri=None):
        """Searches the RDF graph for the annotation ``{namespace_uri}annotation_name``
        for the given symbol and returns whether it has annotation, optionally restricted
        to a specific ``{namespace_uri}annotation_name``

        Specifically, this method searches for a ``annotation_name`` for
        subject symbol
        predicate ``http://biomodels.net/biology-qualifiers/is`` and the object
        specified by ``{namespace_uri}annotation_name``
        for a specific ``{namespace_uri}`` if set, otherwise for any namespace_uri

        Will return a boolean.
        """
        return len(self.get_ontology_terms_by_symbol(symbol, namespace_uri)) != 0

    def get_value(self, symbol):
        """Returns the evaluated value of the given symbol's RHS.
        """
        # Find RHS
        rhs = self.graph.nodes[symbol]['equation'].rhs

        # Evaluate and return
        return float(rhs.evalf())

    def get_initial_value(self, symbol):
        """Returns the initial value of the given symbol

        :param symbol: Sympy Dummy object of required symbol
        :return: float of initial value
        """
        return float(self.dummy_metadata[symbol].initial_value)

    def find_symbols_and_derivatives(self, expression_or_expressions):
        """ Returns a set containing all symbols and derivatives referenced in an expression or list of expressions.

        :param expression_or_expressions: expression or list of expressions to get symbols for
        """
        symbols = set()
        if isinstance(expression_or_expressions, list):
            for expr in expression_or_expressions:
                symbols.update(self.find_symbols_and_derivatives(expr))
        # if this expression is a derivative or a dummy symbol which is not a number placeholder
        elif expression_or_expressions.is_Derivative \
                or (expression_or_expressions.is_Dummy and not self.get_meta_dummy(expression_or_expressions).number):
            symbols.add(expression_or_expressions)
        # otherwise, descend into sub-expressions and collect symbols
        else:
            for arg in expression_or_expressions.args:
                symbols |= self.find_symbols_and_derivatives(arg)

        return symbols

    def find_variable(self, search_dict):
        """Searches for variables by attributes in their MetaDummy record. Pass a dictionary of
        {name: value}, where name is one the attributes in MetaDummy """
        matches = []

        # for each component in this model
        for variable in self.dummy_metadata.values():
            matched = True
            for search_key, search_value in search_dict.items():
                if not (hasattr(variable, search_key) and
                        getattr(variable, search_key) == search_value):
                    matched = False
                    break
            if matched:
                matches.append(variable)
        return matches

    def set_equation(self, lhs, rhs):
        """
        Adds an equation defining the variable named in ``lhs``, or replaces an existing one.

        As with :meth:`add_equation()` the LHS must be either a variable symbol or a derivative, and all numbers and
        variable symbols used in ``lhs`` and ``rhs`` must have been obtained from this model, e.g. via
        :meth:`add_number()`, :meth:`add_variable()`, or :meth:`get_symbol_by_cmeta_id()`.

        :param lhs: An LHS expression (either a symbol or a derivative).
        :param rhs: The new RHS expression for this variable.
        """
        # Get variable symbol named in the lhs
        lhs_symbol = lhs
        if lhs_symbol.is_Derivative:
            lhs_symbol = lhs_symbol.free_symbols.pop()
        assert lhs_symbol.is_Dummy and not self.get_meta_dummy(lhs_symbol).number

        # Check if the variable named in the lhs already has an equation
        i_existing = None
        for i, eq in enumerate(self.equations):
            symbol = eq.lhs.free_symbols.pop() if eq.lhs.is_Derivative else eq.lhs
            if symbol == lhs_symbol:
                i_existing = i
                break

        # Add or replace equation
        if i_existing is None:
            self.equations.append(sympy.Eq(lhs, rhs))
        else:
            self.equations[i_existing] = sympy.Eq(lhs, rhs)

        # Invalidate cache
        self.graph = None

