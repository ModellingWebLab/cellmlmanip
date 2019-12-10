"""Classes to represent a flattened CellML model and metadata about its variables."""
import logging
from io import StringIO

import networkx as nx
import rdflib
import sympy

from cellmlmanip.rdf import create_rdf_node
from cellmlmanip.units import UnitStore


logger = logging.getLogger(__name__)


# Delimiter for variables name in Sympy expressions: <component><delimiter><name>
SYMPY_SYMBOL_DELIMITER = '$'


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
    def __init__(elf, name, cmeta_id=None):

        elf.name = name
        elf.cmeta_id = cmeta_id
        elf.rdf_identity = rdflib.URIRef('#' + cmeta_id) if cmeta_id else None

        # A list of sympy.Eq equation objects
        elf.equations = []

        # A pint UnitStore, mapping unit names to unit objects
        elf.units = UnitStore(model=elf)

        # Maps string variable names to sympy.Dummy objects
        elf._name_to_symbol = dict()

        # Cached nx.DiGraph of this model's equations, with number dummies or with sympy.Number objects
        elf._graph = None
        elf._graph_with_sympy_numbers = None

        # An RDF graph containing further meta data
        elf.rdf = rdflib.Graph()

    def add_unit(elf, name, attributes=None, base_units=False):
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
            elf.units.add_base_unit(name)
        else:
            elf.units.add_custom_unit(name, attributes)

    def add_equation(elf, equation):
        """
        Adds an equation to this model.

        The left-hand side (LHS) of the equation must be either a variable symbol or a derivative.

        All numbers and variable symbols used in the equation must have been obtained from this model, e.g. via
        :meth:`add_number()`, :meth:`add_variable()`, or :meth:`get_symbol_by_ontology_term()`.

        :param equation: A ``sympy.Eq`` object.
        """
        assert isinstance(equation, sympy.Eq), 'The argument `equation` must be a sympy.Eq.'
        elf.equations.append(equation)
        elf._invalidate_cache()

    def add_number(elf, value, units):
        """
        Creates and returns a :class:`NumberDummy` to represent a number with units in sympy expressions.

        :param number: A number (anything convertible to float).
        :param units: A `pint` units representation

        :return: A :class:`NumberDummy` object.
        """

        # Check units
        if not isinstance(units, elf.units.ureg.Unit):
            units = elf.units.get_quantity(units)

        return NumberDummy(value, units)

    def add_variable(elf, name, units, initial_value=None,
                     public_interface=None, private_interface=None, cmeta_id=None):
        """
        Adds a variable to the model and returns a :class:`VariableDummy` to represent it in sympy expressions.

        :param name: A string name.
        :param units: A `pint` units representation.
        :param initial_value: An optional initial value.
        :param public_interface: An optional public interface specifier (only required when parsing CellML).
        :param private_interface: An optional private interface specifier (only required when parsing CellML).
        :param cmeta_id: An optional string specifying a cmeta:id

        :return: A :class:`VariableDummy` object.
        """
        # Check for clashes
        if name in elf._name_to_symbol:
            raise ValueError('Variable %s already exists.' % name)

        # Check units
        if not isinstance(units, elf.units.ureg.Unit):
            units = elf.units.get_quantity(units)

        # Add variable
        elf._name_to_symbol[name] = var = VariableDummy(
            name=name,
            units=units,
            initial_value=initial_value,
            public_interface=public_interface,
            private_interface=private_interface,
            order_added=len(elf._name_to_symbol),
            cmeta_id=cmeta_id,
        )

        # Invalidate cached graphs
        elf._invalidate_cache()

        return var

    def connect_variables(elf, source_name: str, target_name: str):
        """Given the source and target component and variable, create a connection by assigning
        the symbol from the source to the target. If units are not the same, it will add an equation
        to the target component reflecting the relationship. If a symbol has not been assigned to
        the source variable, then return False.

        :param source_name: source variable name
        :param target_name: target variable name
        """
        logger.debug('connect_variables(%s ⟶ %s)', source_name, target_name)

        source = elf._name_to_symbol[source_name]
        target = elf._name_to_symbol[target_name]

        # If the source variable has not been assigned a symbol, we can't make this connection
        if not source.assigned_to:
            logger.info('The source variable has not been assigned to a symbol '
                        '(i.e. expecting a connection): %s ⟶ %s',
                        target.name, source.name)
            return False

        # If target is already assigned this is an error
        if target.assigned_to:
            raise ValueError('Target already assigned to %s before assignment to %s' %
                             (target.assigned_to, source.assigned_to))

        # If source/target variable is in the same unit
        if source.units == target.units:
            # Direct substitution is possible
            target.assigned_to = source.assigned_to
            # everywhere the target variable is used, replace with source variable
            for index, equation in enumerate(elf.equations):
                elf.equations[index] = equation.xreplace({target: source.assigned_to})

        # Otherwise, this connection requires a conversion
        else:
            # Get the scaling factor required to convert source units to target units
            factor = elf.units.convert_to(1 * source.units, target.units).magnitude

            # Dummy to represent this factor in equations, having units for conversion
            factor_dummy = elf.add_number(factor, target.units / source.units)

            # Add an equations making the connection with the required conversion
            elf.equations.append(sympy.Eq(target, source.assigned_to * factor_dummy))

            logger.info('Connection req. unit conversion: %s', elf.equations[-1])

            # The assigned symbol for this variable is itself
            target.assigned_to = target

        logger.debug('Updated target: %s', target)

        # Invalidate cached graphs
        elf._invalidate_cache()

        return True

    def add_rdf(elf, rdf: str):
        """ Takes an RDF string and stores it in the model's RDF graph. """
        elf.rdf.parse(StringIO(rdf))

    def check_left_right_units_equal(elf, equality):
        """
        Checks whether the LHS and RHS in a ``sympy.Eq`` have the same units.
        :param equality: A ``sympy.Eq``.
        """
        lhs_units = elf.units.summarise_units(equality.lhs)
        rhs_units = elf.units.summarise_units(equality.rhs)

        assert elf.units.is_unit_equal(rhs_units, lhs_units), 'Units %s %s != %s %s' % (
            lhs_units, elf.units.ureg.get_base_units(lhs_units),
            rhs_units, elf.units.ureg.get_base_units(rhs_units)
        )

    def get_equations_for(elf, symbols, recurse=True, strip_units=True):
        """Get all equations for a given collection of symbols.

        Results are sorted first by dependencies, then by variable name.

        :param symbols: the symbols to get the equations for.
        :param recurse: indicates whether to recurse the equation graph, or to return only the top level equations.
        :param strip_units: if ``True``, all ``sympy.Dummy`` objects representing number with units will be replaced
            with ordinary sympy number objects.
        """
        # Get graph
        if strip_units:
            graph = elf.graph_with_sympy_numbers
        else:
            graph = elf.graph

        # Get sorted list of symbols
        sorted_symbols = nx.lexicographical_topological_sort(graph, key=str)

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

    def get_derivative_symbols(elf):
        """Returns a list of derivative symbols found in the given model graph.
        The list is ordered by appearance in the cellml document.
        """
        derivative_symbols = [v for v in elf.graph if isinstance(v, sympy.Derivative)]
        return sorted(derivative_symbols, key=lambda state_var: state_var.args[0].order_added)

    def get_state_symbols(elf):
        """Returns a list of state variables found in the given model graph.
        The list is ordered by appearance in the cellml document.
        """
        state_symbols = [v.args[0] for v in elf.get_derivative_symbols()]
        return sorted(state_symbols, key=lambda state_var: state_var.order_added)

    def get_free_variable_symbol(elf):
        """Returns the free variable of the given model graph.
        """
        for v, node in elf.graph.nodes.items():
            if node.get('variable_type', '') == 'free':
                return v

        # This should be unreachable
        raise ValueError('No free variable set in model.')  # pragma: no cover

    def get_rdf_annotations(elf, subject=None, predicate=None, object_=None):
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
        return elf.rdf.triples((subject, predicate, object_))

    def get_rdf_value(elf, subject, predicate):
        """Get the value of an RDF object connected to ``subject`` by ``predicate``.

        :param subject: the object of the triple returned
        :param predicate: the object of the triple returned

        Note: expects exactly one triple to match and the result to be a literal. It's string value is  returned."""
        triples = list(elf.get_rdf_annotations(subject, predicate))
        assert len(triples) == 1
        assert isinstance(triples[0][2], rdflib.Literal)
        value = str(triples[0][2]).strip()  # Could make this cleverer by considering data type if desired
        return value

    def get_symbol_by_cmeta_id(elf, cmeta_id):
        """
        Searches the model and returns the symbol for the variable with the given cmeta id.

        To get symbols from e.g. an oxmeta ontology term, use :meth:`get_symbol_by_ontology_term()`.
        """

        for var in elf._name_to_symbol.values():
            if var.cmeta_id == cmeta_id:
                return var

        raise KeyError('No variable with cmeta id "%s" found.' % str(cmeta_id))

    def get_symbol_by_name(elf, name):
        """ Returns the symbol for the variable with the given ``name``. """
        return elf._name_to_symbol[name]

    def get_symbol_by_ontology_term(elf, namespace_uri, local_name):
        """Searches the RDF graph for a variable annotated with the given
        ``{namespace_uri}local_name`` and returns its symbol.

        Specifically, this method searches for a unique variable annotated with
        predicate ``http://biomodels.net/biology-qualifiers/is`` and the object
        specified by ``{namespace_uri}local_name``.

        Will raise a ``KeyError`` if no variable with the given annotation is
        found, and a ``ValueError`` if more than one variable with the given
        annotation is found.
        """
        symbols = elf.get_symbols_by_rdf(
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

    def get_symbols_by_rdf(elf, predicate, object_=None):
        """Searches the RDF graph for variables annotated with the given predicate and object (e.g. "is oxmeta:time")
        and returns the associated symbols sorted in document order.

        Both ``predicate`` and ``object_`` (if given) must be ``(namespace, local_name)`` tuples or string literals.
        """
        predicate = create_rdf_node(predicate)
        object_ = create_rdf_node(object_)

        # Find symbols
        symbols = []
        for result in elf.rdf.subjects(predicate, object_):
            assert isinstance(result, rdflib.URIRef), 'Non-resource annotated.'

            # Get cmeta id from result uri
            uri = str(result)
            if uri[0] != '#':
                # TODO This should eventually be implemented
                raise NotImplementedError(
                    'Non-local annotations are not supported.')
            symbols.append(elf.get_symbol_by_cmeta_id(uri[1:]))

        return sorted(symbols, key=lambda sym: sym.order_added)

    def get_ontology_terms_by_symbol(elf, symbol, namespace_uri=None):
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
        cmeta_id = elf.graph.nodes[symbol].get('cmeta_id', None)
        if cmeta_id:
            predicate = ('http://biomodels.net/biology-qualifiers/', 'is')
            predicate = create_rdf_node(predicate)
            for delimeter in ('#', '/'):  # Look for terms using either possible namespace delimiter
                subject = rdflib.term.URIRef(delimeter + cmeta_id)
                for object in elf.rdf.objects(subject, predicate):
                    # We are only interested in annotation within the namespace
                    if namespace_uri is None or str(object).startswith(namespace_uri):
                        uri_parts = str(object).split(delimeter)
                        ontology_terms.append(uri_parts[-1])
        return ontology_terms

    def get_units(elf, name):
        """
        Looks up and returns a pint `Unit` object with the given name.
        """
        return elf.units.get_quantity(name)

    @property
    def graph(elf):
        """ A ``networkx.DiGraph`` containing the model equations. """
        # TODO: Set the parameters of the model (parameters rather than use initial values)

        # Return cached graph
        if elf._graph is not None:
            return elf._graph

        # Store symbols, their attributes and their relationships in a directed graph
        graph = nx.DiGraph()

        equation_count = 0

        # Add a node for every variable in the model, and set additional variable meta data
        for equation in elf.equations:
            equation_count += 1

            # Determine LHS.
            lhs = equation.lhs
            if not (lhs.is_Derivative or isinstance(lhs, VariableDummy)):
                raise RuntimeError('DAEs are not supported. All equations must be of form `x = ...` or `dx/dt = ...')

            # Add the lhs symbol of the equation to the graph
            graph.add_node(lhs, equation=equation)

            # Update variable meta data based on the variable's role in the model
            if lhs.is_Derivative:
                # Get the state symbol and update the variable information
                state_symbol = lhs.free_symbols.pop()
                state_symbol.type = 'state'

                # Get the free symbol and update the variable information
                free_symbol = lhs.variables[0]
                free_symbol.type = 'free'
            else:
                lhs.type = None

        # Sanity check: none of the lhs have the same hash
        assert len(graph.nodes) == equation_count

        # Sanity check: all the lhs are unique in meaning (sympy.Dummy: same name != same hash)
        assert len(set([str(x) for x in graph.nodes])) == equation_count

        # Add edges between the nodes
        for equation in elf.equations:
            lhs = equation.lhs

            # for each of the symbols or derivatives on the rhs of the equation
            for rhs in elf.find_symbols_and_derivatives([equation.rhs]):

                if rhs in graph.nodes:
                    # If the symbol maps to a node in the graph just add the dependency edge
                    graph.add_edge(rhs, lhs)
                elif rhs.type in ['state', 'free']:
                    # If the variable is a state or free variable of a derivative
                    graph.add_node(rhs, equation=None, variable_type=rhs.type)
                    graph.add_edge(rhs, lhs)
                else:
                    # this variable is a parameter - add to graph and connect to lhs
                    rhs.type = 'parameter'
                    dummy = elf.add_number(rhs.initial_value, rhs.units)
                    graph.add_node(rhs, equation=sympy.Eq(rhs, dummy), variable_type='parameter')
                    graph.add_edge(rhs, lhs)

            # check that the free  and state variables are defined as a node
            if lhs.is_Derivative:
                state_symbol = lhs.free_symbols.pop()
                free_symbol = lhs.variables[0]
                if free_symbol not in graph.nodes:
                    graph.add_node(free_symbol, equation=None, variable_type=free_symbol.type)
                if state_symbol not in graph.nodes:
                    graph.add_node(state_symbol, equation=None, variable_type=state_symbol.type)

        # Add more meta-data to the graph
        for variable in graph.nodes:
            if not variable.is_Derivative:
                for key in ['cmeta_id', 'name', 'units']:
                    if getattr(variable, key):
                        graph.nodes[variable][key] = getattr(variable, key)
                if variable.type == 'state' and variable.initial_value is not None:
                    graph.nodes[variable]['initial_value'] = sympy.Float(variable.initial_value)
                if variable.type is not None:
                    graph.nodes[variable]['variable_type'] = variable.type

        # Cache graph and return
        elf._graph = graph
        return graph

    @property
    def graph_with_sympy_numbers(elf):
        """
        A ``networkx.DiGraph`` containing the model equations, but with numbers represented as sympy ``Number`` objects
        instead of dummies.
        """
        if elf._graph_with_sympy_numbers is not None:
            return elf._graph_with_sympy_numbers

        # Get a clone of the graph
        graph = elf.graph.copy()

        # Replace dummies with Float objects
        for node in graph.nodes:
            equation = graph.nodes[node]['equation']
            if equation is None:
                continue

            # Get all the dummy symbols on the RHS
            dummies = equation.rhs.atoms(sympy.Dummy)

            # Get any dummy symbols which are placeholders for numbers
            subs_dict = {d: float(d) for d in dummies if isinstance(d, NumberDummy)}

            # And replace the equation with one with the rhs subbed with sympy.Number objects
            if subs_dict:
                # Update rhs
                rhs = equation.rhs.subs(subs_dict)

                # Check if simplification removed dependencies on other variables, and if so remove the corresponding
                # edges.
                refs = elf.find_symbols_and_derivatives([rhs])
                edges = list(graph.in_edges(equation.lhs))
                for edge in edges:
                    ref = edge[0]
                    if ref not in refs:
                        graph.remove_edge(ref, equation.lhs)

                # Replace equation
                graph.nodes[node]['equation'] = sympy.Eq(equation.lhs, rhs)

        # Cache graph and return
        elf._graph_with_sympy_numbers = graph
        return graph

    def has_ontology_annotation(elf, symbol, namespace_uri=None):
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
        return len(elf.get_ontology_terms_by_symbol(symbol, namespace_uri)) != 0

    def _invalidate_cache(elf):
        """ Removes cached graphs: should be called after manipulating variables or equations. """
        elf._graph = None
        elf._graph_with_sympy_numbers = None

    def get_value(elf, symbol):
        """ Returns the evaluated value of the given symbol's RHS. """
        return float(elf.graph.nodes[symbol]['equation'].rhs.evalf())

    def get_initial_value(elf, symbol):
        """
        Returns the initial value of the given symbol.

        :param symbol: Sympy Dummy object of required symbol
        :return: float of initial value
        """
        return symbol.initial_value

    def find_symbols_and_derivatives(elf, expression):
        """ Returns a set containing all symbols and derivatives referenced in a list of expressions.

        :param expression: a list of expressions to get symbols for.
        :return: a set of symbols and derivative objects.
        """
        symbols = set()
        for expr in expression:
            if expr.is_Derivative or isinstance(expr, VariableDummy):
                symbols.add(expr)
            else:
                symbols |= elf.find_symbols_and_derivatives(expr.args)
        return symbols

    def remove_equation(elf, equation):
        """
        Removes an equation from the model.

        :param equation: The equation to remove.
        """
        try:
            elf.equations.remove(equation)
        except ValueError:
            raise KeyError('Equation not found in model ' + str(equation))

        # Invalidate cached equation graphs
        elf._invalidate_cache()

    def variables(elf):
        """ Returns an iterator over this model's variable symbols. """
        return elf._name_to_symbol.values()


class NumberDummy(sympy.Dummy):
    """
    Used to represent a number with a unit, inside a Sympy expression.

    Unlike sympy expressions, this number type will never be removed in simplify operations etc.

    Number dummies should never be created directly, but always via :meth:`Model.add_number()`.
    """
    # Sympy annoyingly overwrites __new__
    def __new__(cls, value, *args, **kwargs):
        return super().__new__(cls, str(value))

    def __init__(elf, value, units):
        elf.value = float(value)
        elf.units = units

    def __float__(elf):
        return elf.value

    def __str__(elf):
        return str(elf.value)


class VariableDummy(sympy.Dummy):
    """
    Used to represent a variable (with meta data) in a Sympy expression.

    Variable dummies should never be created directly, but always via :meth:`Model.add_variable()`.

    For the constructor arguments, see :meth:`Model.add_variable()`.
    """
    # Sympy annoyingly overwrites __new__
    def __new__(cls, name, *args, **kwargs):
        return super().__new__(cls, name)

    def __init__(elf,
                 name,
                 units,
                 initial_value=None,
                 public_interface=None,
                 private_interface=None,
                 order_added=None,
                 cmeta_id=None):

        elf.name = name
        elf.units = units
        elf.initial_value = None if initial_value is None else float(initial_value)

        # Interface properties, only used during parsing
        elf.public_interface = public_interface
        elf.private_interface = private_interface

        # Variables are either 'source' variables, or receive their value from
        # a variable that they're connected to (using CellML connections).
        # The ``assigned_to`` property is used to indicate where this object
        # receives its value.
        elf.assigned_to = None
        if not (private_interface == 'in' or public_interface == 'in'):
            elf.assigned_to = elf

        # Optional order added, used for sorting sometimes.
        elf.order_added = order_added

        # Optional cmeta id
        elf.cmeta_id = cmeta_id

        # This variable's type
        # TODO: Define allowed types via enum
        elf.type = None

    def __str__(elf):
        return elf.name

