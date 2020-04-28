"""Classes to represent a flattened CellML model and metadata about its variables."""
import logging
from enum import Enum
from io import StringIO

import networkx as nx
import rdflib
import sympy

from cellmlmanip.rdf import create_rdf_node
from cellmlmanip.units import UnitStore


logger = logging.getLogger(__name__)

# Delimiter for variables name in Sympy expressions: <component><delimiter><name>
SYMPY_SYMBOL_DELIMITER = '$'

# Float precision to use when creating sympy.Float objects
FLOAT_PRECISION = 17


class DataDirectionFlow(Enum):
    """ Direction of data flow for converting units"""
    INPUT = 1
    OUTPUT = 2


class Model(object):
    """
    A componentless representation of a CellML model, containing a list of equations, units, and RDF metadata about
    variables used in those equations.

    The main parts of a Model are 1. a list of sympy equation objects; 2. a collection of named units; and 3. an RDF
    graph that stores further meta data about the model.

    Equations are stored as ``Sympy.Eq`` objects, but with the caveat that all variables and numbers must be specified
    using the ``Sympy.Dummy`` objects returned by :meth:`add_variable()` and :meth:`add_number()`.

    Units are handled using the ``units`` property of a model, which is an instance of ``cellmlmanip.units.UnitStore``.

    RDF meta data can be attached to either the model itself or to model variables, via the CellML ``cmeta:id``
    attribute. Cmeta ids set on any other parts of CellML models are ignored.

    Cellmlmanip does not support algebraic models: the left-hand side every equation in the model must be a variable or
    a derivative.

    :param name: the name of the model e.g. from ``<model name="">``.
    :param cmeta_id: An optional cmeta id, e.g. from ``<model cmeta:id="">``.
    :param unit_store: Optional :class:`cellmlmanip.units.UnitStore` instance; if given the model will share the
        underlying registry so that conversions between model units and those from the provided store work.
    """

    def __init__(self, name, cmeta_id=None, unit_store=None):

        self.name = name
        self._cmeta_id = cmeta_id
        self.rdf_identity = create_rdf_node('#' + cmeta_id) if cmeta_id else None

        # A list of sympy.Eq equation objects
        self.equations = []

        # A UnitStore object
        if unit_store:
            self.units = UnitStore(unit_store)
        else:
            self.units = UnitStore()

        # Maps string variable names to VariableDummy objects
        self._name_to_variable = {}

        # Maps cmeta ids to VariableDummy objects
        self._cmeta_id_to_variable = {}

        # Cached nx.DiGraph of this model's equations, with number dummies or with sympy.Number objects
        self._graph = None
        self._graph_with_sympy_numbers = None

        # An RDF graph containing further meta data
        self.rdf = rdflib.Graph()

        # Map from VariableDummy to defining equation, where the variable is defined by a simple equation
        self._var_definition_map = {}

        # Map from VariableDummy to defining equation, where the variable is defined by an ODE
        self._ode_definition_map = {}

    def add_cmeta_id(self, variable):
        """
        Adds a (unique) cmeta id to the given variable.

        If the variable already has a cmeta id no action is performed.

        :param variable: A :class:`VariableDummy`.
        """
        if variable._cmeta_id is not None:
            return

        # Create new cmeta id
        cmeta_id = str(variable)
        while self.has_cmeta_id(cmeta_id):
            cmeta_id += '_'

        # Add to variable and store in mapping
        variable._set_cmeta_id(cmeta_id)
        self._cmeta_id_to_variable[cmeta_id] = variable

    def add_equation(self, equation, check_duplicates=True):
        """
        Adds an equation to this model.

        The left-hand side (LHS) of the equation must be either a variable (as a :class:`VariableDummy`) or a derivative
        (as a ``sympy.Derivative``).

        All numbers and variables used in the equation must have been obtained from this model, e.g. via
        :meth:`add_number()`, :meth:`add_variable()`, or :meth:`get_variable_by_ontology_term()`.

        :param equation: A ``sympy.Eq`` object.
        :param check_duplicates: whether to check that the equation's LHS is not already defined
        """
        assert isinstance(equation, sympy.Eq), 'The argument `equation` must be a sympy.Eq.'
        self.equations.append(equation)
        lhs = equation.lhs
        if lhs.is_Derivative:
            state_var = lhs.free_symbols.pop()
            if check_duplicates:
                self._check_duplicate_definitions(state_var, equation)
            self._ode_definition_map[state_var] = equation
        elif isinstance(lhs, VariableDummy):
            if check_duplicates:
                self._check_duplicate_definitions(lhs, equation)
            self._var_definition_map[lhs] = equation
        else:
            raise ValueError('Equation LHS should be a derivative or variable, not {}'.format(lhs))
        self._invalidate_cache()

    def _check_duplicate_definitions(self, var, equation):
        """
        Assert that a variable doesn't have an existing definition.

        :param var: the VariableDummy, either a state var or normal var
        :param equation: the new definition being added
        """
        if var in self._ode_definition_map:
            raise ValueError('The variable {} is defined twice ({} and {})'.format(
                var, equation, self._ode_definition_map[var]))
        if var in self._var_definition_map:
            raise ValueError('The variable {} is defined twice ({} and {})'.format(
                var, equation, self._var_definition_map[var]))

    def add_number(self, value, units):
        """
        Creates and returns a :class:`NumberDummy` to represent a number with units in sympy expressions.

        :param number: A number (anything convertible to float).
        :param units: A string unit name or a ``Unit`` object.

        :return: A :class:`NumberDummy` object.
        """

        # Check units
        if not isinstance(units, self.units.Unit):
            units = self.units.get_unit(units)

        return NumberDummy(value, units)

    def add_variable(self, name, units, initial_value=None, public_interface=None, private_interface=None,
                     cmeta_id=None):
        """
        Adds a variable to the model and returns a :class:`VariableDummy` to represent it in sympy expressions.

        :param name: A string name.
        :param units: A string unit name or a ``Unit`` object.
        :param initial_value: An optional initial value.
        :param public_interface: An optional public interface specifier (only required when parsing CellML).
        :param private_interface: An optional private interface specifier (only required when parsing CellML).
        :param cmeta_id: An optional string specifying a cmeta:id
        :raises ValueError: If a variable with that name already exists, or the given cmeta_id is already taken.
        :return: A :class:`VariableDummy` object.
        """
        # Check for clashes
        if name in self._name_to_variable:
            raise ValueError('Variable %s already exists.' % name)

        # Check uniqueness of cmeta id
        if cmeta_id is not None and self.has_cmeta_id(cmeta_id):
            raise ValueError('The cmeta id "%s" is already in use.' % cmeta_id)

        # Check units
        if not isinstance(units, self.units.Unit):
            units = self.units.get_unit(units)

        # Add variable
        self._name_to_variable[name] = var = VariableDummy(
            name=name,
            units=units,
            model=self,
            initial_value=initial_value,
            public_interface=public_interface,
            private_interface=private_interface,
            order_added=len(self._name_to_variable),
            cmeta_id=cmeta_id,
        )

        # Add cmeta id to var mapping
        if cmeta_id is not None:
            self._cmeta_id_to_variable[cmeta_id] = var

        # Invalidate cached graphs
        self._invalidate_cache()

        return var

    def transform_constants(self):
        """
        Called by CellML parser to standardise handling of 'constants'.

        Once this has been called, the only variables with an initial_value attribute will be state variables,
        and the initial value will do what it implies - hold the value the state variable should take at t=0.

        Non state variables with an initial value are actually just constants. For consistent processing later on we add
        equations defining them, and remove the initial_value attribute.
        """
        for var in self._name_to_variable.values():
            if var in self._ode_definition_map:
                assert var.initial_value is not None, 'State variable {} has no initial_value set'.format(var)
            elif var.initial_value is not None:
                value = self.add_number(var.initial_value, var.units)
                self.add_equation(sympy.Eq(var, value))
                var.initial_value = None

    def connect_variables(self, source_name: str, target_name: str):
        """Tells this model that the variable indicated by ``target_name`` should get its value from the variable
        indicated by ``source_name``.

        If the variables' units are not the same, an equation will be added to the target component.

        This method will only work if the source variable has been assigned a value, either through an equation or by
        connecting it (directly or indirectly) to a variable with an equation. If this is not yet the case, ``False`` is
        returned. If successful the method returns ``True``.

        :param source_name: The source variable name
        :param target_name: The target variable name
        """
        logger.debug('connect_variables(%s ⟶ %s)', source_name, target_name)

        source = self._name_to_variable[source_name]
        target = self._name_to_variable[target_name]

        # If the source variable has not been assigned a value, we can't make this connection
        if not source.assigned_to:
            logger.debug(
                'Cannot connect {} to {} at this time: The source variable has not been assigned a value '.format(
                    target.name, source.name))
            return False

        # If target is already assigned this is an error
        if target.assigned_to:
            raise ValueError('Target already assigned to %s before assignment to %s' %
                             (target.assigned_to, source.assigned_to))

        if target in self._var_definition_map:
            raise ValueError('Multiple definitions for {} ({} and {})'.format(
                target, source, self._var_definition_map[target]))

        # If source/target variable is in the same unit
        if source.units == target.units:
            # Direct substitution is possible
            target.assigned_to = source.assigned_to
            # everywhere the target variable is used, replace with source variable,
            # updating the definition maps accordingly
            for index, equation in enumerate(self.equations):
                self.equations[index] = new_eq = equation.xreplace({target: source.assigned_to})
                if equation.lhs.is_Derivative:
                    state_var = equation.lhs.free_symbols.pop()
                    assert state_var in self._ode_definition_map
                    self._ode_definition_map[state_var] = new_eq
                else:
                    assert equation.lhs in self._var_definition_map
                    self._var_definition_map[equation.lhs] = new_eq
            if target in self._ode_definition_map:
                self._ode_definition_map[source.assigned_to] = self._ode_definition_map[target]
                del self._ode_definition_map[target]
            # Migrate the cmeta:id too so annotations work as expected
            # Really modellers should annotate the source variable, but they don't always!
            if target.cmeta_id is not None:
                # Note that this method considers source as the variable with the cmeta:id already.
                self.transfer_cmeta_id(source=target, target=source)

        # Otherwise, this connection requires a conversion
        else:
            # Get the scaling factor required to convert source units to target units
            factor = self.units.convert(1 * source.units, target.units).magnitude

            # Dummy to represent this factor in equations, having units for conversion
            factor_dummy = self.add_number(factor, target.units / source.units)

            # Add an equations making the connection with the required conversion
            self.add_equation(sympy.Eq(target, source.assigned_to * factor_dummy))

            logger.info('Connection req. unit conversion: %s', self.equations[-1])

            # The assigned variable for this variable is itself
            target.assigned_to = target

        logger.debug('Updated target: %s', target)

        # Invalidate cached graphs
        self._invalidate_cache()

        return True

    def add_rdf(self, rdf: str):
        """ Takes an RDF string and stores it in the model's RDF graph. """
        self.rdf.parse(StringIO(rdf))

    def check_left_right_units_equal(self, equality):
        """
        Checks whether the LHS and RHS in a ``sympy.Eq`` have the same units.
        :param equality: A ``sympy.Eq``.
        """
        lhs_units = self.units.evaluate_units(equality.lhs)
        rhs_units = self.units.evaluate_units(equality.rhs)
        assert self.units.is_equivalent(rhs_units, lhs_units), 'Units %s %s != %s %s' % (
            self.units.format(lhs_units), self.units.format(lhs_units, True),
            self.units.format(rhs_units), self.units.format(rhs_units, True),
        )

    def get_equations_for(self, variables, recurse=True, strip_units=True):
        """Get all equations for a given collection of variables.

        Results are sorted first by dependencies, then by variable name.

        :param variables: The variables to get the equations for (as :class:`VariableDummy` objects).
        :param recurse: Indicates whether to recurse the equation graph, or to return only the top level equations.
        :param strip_units: If ``True``, all ``sympy.Dummy`` objects representing number with units will be replaced
            with ordinary sympy number objects.
        """
        # Get graph
        if strip_units:
            graph = self.graph_with_sympy_numbers
        else:
            graph = self.graph

        # Get sorted list of variables
        sorted_variables = nx.lexicographical_topological_sort(graph, key=str)

        # Create set of variables for which we require equations
        required_variables = set()
        for output in variables:
            required_variables.add(output)
            if recurse:
                required_variables.update(nx.ancestors(graph, output))
            else:
                required_variables.update(graph.pred[output])

        eqs = []
        for variable in sorted_variables:
            # Ignore variables we don't need
            if variable not in required_variables:
                continue

            # Get equation
            eq = graph.nodes[variable]['equation']

            # Skip variables that are not set with an equation
            if eq is None:
                continue

            eqs.append(eq)

        return eqs

    def get_derivatives(self):
        """Returns a list of ``sympy.Derivative`` objects found as LHS in the given model graph.
        The list is ordered by appearance in the cellml document.
        """
        derivatives = [v for v in self.graph if isinstance(v, sympy.Derivative)]
        return sorted(derivatives, key=lambda deriv: deriv.args[0].order_added)

    def get_derived_quantities(self):
        """Returns a list of derived quantities found in the given model graph.
        A derived quantity is any variable that is not a state variable or parameter/constant.
        """
        derived_quantities = [v for v, node in self.graph.nodes.items()
                              if not isinstance(v, sympy.Derivative)
                              and node.get('variable_type', '') not in ('state', 'free', 'parameter')]
        return sorted(derived_quantities, key=lambda var: var.order_added)

    def get_display_name(self, var, ontology=None):
        """Return a display name for the given variable.

        Looks for an annotation in the ontology first (or the local name from any annotation if no ontology is
        specified), then cmeta:id if present, or the variable's name attribute if not.
        
        Dollar symbols in the name are replaced by a double underscore.
        :param var: the variable for which to get the display name.
        :param ontology: the namespace prefix for the ontology used (defaults to not considering annotations).

        :return: the display name for the variable according to the algorithm above
        """
        if self.has_ontology_annotation(var, ontology):
            display_name = self.get_ontology_terms_by_variable(var, ontology)[-1]
        else:
            display_name = var.cmeta_id if var.cmeta_id else var.name
        return display_name.replace('$', '__')

    def get_state_variables(self):
        """
        Returns a list of state variables found in the given model graph (ordered by appearance in the CellML document).
        """
        states = list(self._ode_definition_map.keys())
        return sorted(states, key=lambda state_var: state_var.order_added)

    def get_free_variable(self):
        """Returns the free variable in this model (if any)."""
        for v, node in self.graph.nodes.items():
            if node.get('variable_type', '') == 'free':
                return v

        raise ValueError('No free variable set in model.')  # pragma: no cover

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

    def get_variable_by_cmeta_id(self, cmeta_id):
        """
        Searches the model and returns the variable with the given cmeta id.

        To get variables from e.g. an oxmeta ontology term, use :meth:`get_variable_by_ontology_term()`.

        :param cmeta_id: Either a string id or :class:`rdflib.URIRef` instance.
        :returns: A :class:`VariableDummy` object
        """

        # Get cmeta_id from URIRef
        if isinstance(cmeta_id, rdflib.term.Node):
            assert isinstance(cmeta_id, rdflib.URIRef), 'Non-resource {} annotated.'.format(cmeta_id)
            cmeta_id = str(cmeta_id)
            if cmeta_id[0] != '#':
                # TODO This should eventually be implemented?
                raise NotImplementedError(
                    'Non-local annotations are not supported.')
            cmeta_id = cmeta_id[1:]

        # Get variable by cmeta id string
        assert isinstance(cmeta_id, str)
        try:
            return self._cmeta_id_to_variable[cmeta_id]
        except KeyError:
            raise KeyError('No variable with cmeta id "%s" found.' % str(cmeta_id))

    def get_variable_by_name(self, name):
        """ Returns the variable with the given ``name``. """
        return self._name_to_variable[name]

    def get_variable_by_ontology_term(self, term):
        """Searches the RDF graph for a variable annotated with the given ``term`` and returns it.

        Specifically, this method searches for a unique variable annotated with
        predicate ``http://biomodels.net/biology-qualifiers/is`` and the object
        specified by ``term``.

        Will raise a ``KeyError`` if no variable with the given annotation is
        found, and a ``ValueError`` if more than one variable with the given
        annotation is found.

        :param term: anything suitable as an input to :meth:`create_rdf_node`; typically either an RDF
            node already, or a tuple ``(namespace_uri, local_name)``.
        """
        variables = self.get_variables_by_rdf(('http://biomodels.net/biology-qualifiers/', 'is'), term)
        if len(variables) == 1:
            return variables[0]
        elif not variables:
            raise KeyError('No variable annotated with {} found.'.format(term))
        else:
            raise ValueError('Multiple variables annotated with {}'.format(term))

    def get_variables_by_rdf(self, predicate, object_=None):
        """Searches the RDF graph for variables annotated with the given predicate and object (e.g. ``is oxmeta:time``)
        and returns the associated variables sorted in document order.

        Both ``predicate`` and ``object_`` (if given) must be ``(namespace, local_name)`` tuples or string literals.
        """
        predicate = create_rdf_node(predicate)
        object_ = create_rdf_node(object_)

        # Find variables, sort and return
        variables = [self.get_variable_by_cmeta_id(result) for result in self.rdf.subjects(predicate, object_)]
        return sorted(variables, key=lambda sym: sym.order_added)

    def get_ontology_terms_by_variable(self, variable, namespace_uri=None):
        """
        Returns all ontology terms linked to the variable ``variable`` via the
        ``http://biomodels.net/biology-qualifiers/is`` predicate.

        :param variable: The variable to search for (as a :class:`VariableDummy` object).
        :param namespace_uri: An optional namespace URI. If given, only terms within the given namespace will be
            returned.
        :returns: A list of term names.
        """
        ontology_terms = []
        if variable.rdf_identity:
            predicate = create_rdf_node(('http://biomodels.net/biology-qualifiers/', 'is'))
            for object in self.rdf.objects(variable.rdf_identity, predicate):
                # Filter by namespace
                if namespace_uri is None or str(object).startswith(namespace_uri):
                    uri_parts = str(object).split('#')
                    ontology_terms.append(uri_parts[-1])
        return ontology_terms

    def get_definition(self, variable):
        """Get the equation (if any) defining the given variable.

        :param variable: The variable to look up (as a class:`VariableDummy`. If this appears as the LHS of a straight
            assignment, or the state variable in an ODE, the corresponding equation will be returned.
        :returns: A Sympy equation, or ``None`` if the variable is not defined by an equation.
        """
        defn = self._ode_definition_map.get(variable)
        if defn is None:
            defn = self._var_definition_map.get(variable)
        return defn

    def get_value(self, variable):
        """ Returns the evaluated value of the given variable's RHS. """
        return float(self.graph.nodes[variable]['equation'].rhs.evalf())

    @property
    def graph(self):
        """ A ``networkx.DiGraph`` containing the model equations. """
        # TODO: Set the parameters of the model (parameters rather than use initial values)

        # Return cached graph
        if self._graph is not None:
            return self._graph

        # Store variables, their attributes and their relationships in a directed graph
        graph = nx.DiGraph()

        equation_count = 0

        # Add a node for every variable in the model, and set variable types
        for equation in self.equations:
            equation_count += 1

            # Add the lhs of the equation to the graph
            lhs = equation.lhs
            graph.add_node(lhs, equation=equation)

            # Update variable meta data based on the variable's role in the model
            if lhs.is_Derivative:
                # Get the state symbol and update the variable information
                state_symbol = lhs.free_symbols.pop()
                state_symbol.type = 'state'

                # Get the free symbol and update the variable information
                free_symbol = lhs.variables[0]
                free_symbol.type = 'free'
            elif isinstance(equation.rhs, NumberDummy):
                lhs.type = 'parameter'
            else:
                lhs.type = 'computed'

        # Sanity check: none of the lhs have the same hash
        assert len(graph.nodes) == equation_count

        # Sanity check: all the lhs are unique in meaning (sympy.Dummy: same name != same hash)
        assert len(set([str(x) for x in graph.nodes])) == equation_count

        # Add edges between the nodes
        for equation in self.equations:
            lhs = equation.lhs

            # for each of the symbols or derivatives on the rhs of the equation
            for rhs in self.find_variables_and_derivatives([equation.rhs]):

                if rhs in graph.nodes:
                    # If the symbol maps to a node in the graph just add the dependency edge
                    graph.add_edge(rhs, lhs)
                elif rhs.type in ['state', 'free']:
                    # If the variable is a state or free variable of a derivative
                    graph.add_node(rhs, equation=None, variable_type=rhs.type)
                    graph.add_edge(rhs, lhs)
                else:
                    assert False, 'Unexpected variable {} on RHS'.format(rhs)  # pragma: no cover

            # check that the free  and state variables are defined as a node
            if lhs.is_Derivative:
                state_symbol = lhs.free_symbols.pop()
                free_symbol = lhs.variables[0]
                if free_symbol not in graph.nodes:
                    graph.add_node(free_symbol, equation=None, variable_type=free_symbol.type)
                if state_symbol not in graph.nodes:
                    graph.add_node(state_symbol, equation=None, variable_type=state_symbol.type)

        # Store variable type in the graph too
        for variable in graph.nodes:
            if not variable.is_Derivative:
                assert variable.type is not None
                graph.nodes[variable]['variable_type'] = variable.type

        # Cache graph and return
        self._graph = graph
        return graph

    @property
    def graph_with_sympy_numbers(self):
        """
        A ``networkx.DiGraph`` containing the model equations, but with numbers represented as sympy ``Number`` objects
        instead of dummies.
        """
        if self._graph_with_sympy_numbers is not None:
            return self._graph_with_sympy_numbers

        # Get a clone of the graph
        graph = self.graph.copy()

        # Replace dummies with Float objects
        for node in graph.nodes:
            equation = graph.nodes[node]['equation']
            if equation is None:
                continue

            # Get all the dummy numbers on the RHS, and prepare substitution map
            dummies = equation.rhs.atoms(NumberDummy)
            subs_dict = {d: sympy.Float(d.value, FLOAT_PRECISION) for d in dummies}

            # And replace the equation with one with the rhs subbed with sympy.Number objects
            if subs_dict:
                # Update rhs
                rhs = equation.rhs.xreplace(subs_dict)

                # Check if simplification removed dependencies on other variables, and if so remove the corresponding
                # edges.
                refs = self.find_variables_and_derivatives([rhs])
                edges = list(graph.in_edges(equation.lhs))
                for edge in edges:
                    ref = edge[0]
                    if ref not in refs:
                        graph.remove_edge(ref, equation.lhs)

                # Replace equation
                graph.nodes[node]['equation'] = sympy.Eq(equation.lhs, rhs)

        # Cache graph and return
        self._graph_with_sympy_numbers = graph
        return graph

    def has_cmeta_id(self, cmeta_id):
        """
        Returns ``True`` only if the given ``cmeta_id`` exists in this model.

        Note that only cmeta ids on variables or the model itself are checked and supported.
        """
        # Check if it's the model id
        if cmeta_id == self._cmeta_id and cmeta_id is not None:
            return True

        # Check if it's a variable id. All other cmeta_ids in the original CellML are ignored.
        return cmeta_id in self._cmeta_id_to_variable

    def has_ontology_annotation(self, variable, namespace_uri=None):
        """Checks that there is at least one result for
        :meth:`Model.get_ontology_terms_by_variable(variable, namespace_uri)`.
        """
        return len(self.get_ontology_terms_by_variable(variable, namespace_uri)) != 0

    def _invalidate_cache(self):
        """ Removes cached graphs: should be called after manipulating variables or equations. """
        self._graph = None
        self._graph_with_sympy_numbers = None

    def is_state(self, variable):
        """ Checks if the given ``variable`` is a state variable (i.e. if it's defined by an ODE). """
        return variable in self._ode_definition_map

    def is_constant(self, variable):
        """Determine whether the given ``variable`` is a constant.

        This is calculated by looking at the RHS of the defining equation and checking it has no
        variable references.
        """
        defn = self._var_definition_map.get(variable)
        return defn is not None and len(defn.rhs.atoms(VariableDummy)) == 0

    def find_variables_and_derivatives(self, expressions):
        """ Returns a set containing all variables and derivatives referenced in a list of expressions.

        Note that we can't just use ``.atoms(VariableDummy, sympy.Derivative)`` for this, because it
        will return the state and free variables from inside derivatives, which is not what we want.

        :param expressions: an iterable of expressions to get variables for.
        :return: a set of variables and derivatives, as :class:`VariableDummy` and ``sympy.Derivative`` respectively.
        """
        variables = set()
        for expr in expressions:
            if expr.is_Derivative or isinstance(expr, VariableDummy):
                variables.add(expr)
            else:
                variables |= self.find_variables_and_derivatives(expr.args)
        return variables

    def remove_equation(self, equation):
        """
        Removes an equation from the model.

        :param equation: The equation to remove.
        """
        try:
            self.equations.remove(equation)
        except ValueError:
            raise KeyError('Equation not found in model ' + str(equation))

        # Update dependency maps
        lhs = equation.lhs
        if lhs.is_Derivative:
            del self._ode_definition_map[lhs.free_symbols.pop()]
        else:
            del self._var_definition_map[lhs]

        # Invalidate cached equation graphs
        self._invalidate_cache()

    def remove_variable(self, variable):
        """Remove a variable and its defining equation from the model.

        This will remove the equation either that defines ``variable`` directly, or if it is a state variable the
        corresponding ODE. All annotations about this variable are also removed from the RDF graph.

        :param VariableDummy variable: the variable to remove
        """
        # Remove defining equation
        defn = self.get_definition(variable)
        if defn is not None:
            self.remove_equation(defn)

        # Remove any annotations
        if variable.rdf_identity:
            for triple in self.rdf.triples((variable.rdf_identity, None, None)):
                self.rdf.remove(triple)

        # Remove references to variable and invalidate cache
        del self._name_to_variable[variable.name]
        if variable._cmeta_id is not None:
            del self._cmeta_id_to_variable[variable._cmeta_id]
        variable._model = None  # Just in case!
        self._invalidate_cache()

    def transfer_cmeta_id(self, source, target):
        """
        Removes the ``cmeta_id`` from the variable ``source`` and adds it to ``target``.

        Raises a ``ValueError`` if ``source`` doesn't have a cmeta id, or if ``target`` already has a cmeta id.
        """
        if source._cmeta_id is None:
            raise ValueError('Cannot transfer cmeta id: source variable has no cmeta id.')
        if target._cmeta_id is not None:
            raise ValueError('Cannot transfer cmeta id: target variable already has a cmeta id.')

        # Transfer id
        target._set_cmeta_id(source._cmeta_id)
        source._set_cmeta_id(None)

        # Update mapping
        self._cmeta_id_to_variable[target._cmeta_id] = target

    def variables(self):
        """ Returns an iterator over this model's variables. """
        return self._name_to_variable.values()

    def convert_variable(self, original_variable, units, direction):
        """
        Add a new linked version of the given variable in the desired units.

        If the variable already has the requested units, no changes are made and the original variable is returned.
        Otherwise ``direction`` specifies how information flows between the new variable and the original, and
        hence what new equation(s) are added to the model to perform the conversion.
        If ``INPUT`` then the original variable takes its value from the newly added variable;
        if ``OUTPUT`` then the opposite happens.

        Any ``cmeta:id`` attribute on the original variable is moved to the new one, so ontology annotations will refer
        to the new variable.

        Similarly if the direction is ``INPUT`` then any initial value will be moved to the new variable
        (and converted appropriately).

        For example::

            Original model
                var{time} time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};

                ode(sv1, time) = 1{mV_per_ms};

        convert_variable(time, second, DataDirectionFlow.INPUT)

        becomes
                var time: ms;
                var{time} time_converted: s;
                var{sv11} sv1: mV {init: 2};
                var sv1_orig_deriv mV_per_ms

                time = 1000 * time_converted;
                sv1_orig_deriv = 1{mV_per_ms}
                ode(sv1, time_converted) = 1000 * sv1_orig_deriv


        convert_variable(time, second, DataDirectionFlow.OUTPUT)

            creates model
                var{time} time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};
                var{time} time_converted: s

                ode(sv1, time) = 1{mV_per_ms};
                time_converted = 0.001 * time

        :param original_variable: the VariableDummy object representing the variable in the model to be converted
        :param units: a Pint unit object representing the units to convert variable to (note if variable is already
                      in these units, model remains unchanged and the original variable is returned
        :param direction: either DataDirectionFlow.INPUT; the variable to be changed is an input and all affected
                          equations will be adjusted
                          or DataDirectionFlow.OUTPUT; the variable to be changed is an output, equations
                          are unaffected apart from converting the actual output
        :return: new variable with desired units, or original unchanged if conversion was not necessary
        :raises DimensionalityError: if the unit conversion is impossible
        """
        # assertion errors will be thrown here if arguments are incorrect type
        self._check_arguments_for_convert_variables(original_variable, units, direction)

        original_units = original_variable.units
        # no conversion necessary
        if 1 * original_units == 1 * units:
            return original_variable

        # conversion_factor for old units to new
        # throws DimensionalityError if unit conversion is not possible
        cf = self.units.get_conversion_factor(from_unit=original_units, to_unit=units)

        state_symbols = self.get_state_variables()
        free_symbol = self.get_free_variable()
        # create new variable and relevant equations
        new_variable = self._convert_variable_instance(original_variable, cf, units, direction)

        # if is output do not need to do additional changes for state/free symbols
        if direction == DataDirectionFlow.OUTPUT:
            return new_variable

        new_derivatives = []
        # if state variable
        if original_variable in state_symbols:
            new_derivatives.append(self._convert_state_variable_deriv(original_variable, new_variable, cf))

        # if free variable
        if original_variable == free_symbol:
            # for each derivative wrt to free variable add necessary variables/equations
            for ode in [eq for v, eq in sorted(self._ode_definition_map.items(),
                                               key=lambda v_eq: v_eq[0].order_added)]:
                if ode.args[0].args[1].args[0] == original_variable:
                    new_derivatives.append(self._convert_free_variable_deriv(ode, new_variable, cf))

        # replace any instances of derivative of rhs of other eqns with new derivative variable
        for new_derivative in new_derivatives:
            self._replace_derivatives(new_derivative)

        self._invalidate_cache()

        return new_variable

    def _check_arguments_for_convert_variables(self, variable, units, direction):
        """
        Checks the arguments of the convert_variable function.
        :param variable: variable must be a VariableDummy object present in the model
        :param units: units must be a pint Unit object in this model
        :param direction: must be part of DataDirectionFlow enum
       """
        # variable should be a VariableDummy
        assert isinstance(variable, VariableDummy)

        # variable must be in model
        assert variable.name in self._name_to_variable

        # units should be a pint Unit object in the registry for this model
        assert isinstance(units, self.units.Unit)

        # direction should be part of enum
        assert isinstance(direction, DataDirectionFlow)

    def _replace_derivatives(self, new_derivative):
        """
        Function to replace an instance of a derivative that occurs on the RHS of any equation
        :param new_derivative: new variable representing the derivative
        """
        for equation in self.equations.copy():
            for argument in equation.rhs.atoms(sympy.Derivative):
                if new_derivative['expression'] == argument:
                    # add new equation
                    new_eqn = equation.xreplace({new_derivative['expression']: new_derivative['variable']})
                    self.remove_equation(equation)
                    self.add_equation(new_eqn)
                    break

    def _create_new_deriv_variable_and_equation(self, eqn, derivative_variable):
        """
        Create a new variable and equation for the derivative.
        :param eqn: the original derivative eqn
        :param derivative_variable: the dependent variable
        :return: new variable for the derivative
        """
        # 1. create a new variable
        deriv_name = self.get_unique_name(derivative_variable.name + '_orig_deriv')
        deriv_units = self.units.evaluate_units(eqn.args[0])
        new_deriv_variable = self.add_variable(name=deriv_name, units=deriv_units)

        # 2. create new equation and remove original
        expression = sympy.Eq(new_deriv_variable, eqn.args[1])
        self.remove_equation(eqn)
        self.add_equation(expression)

        return new_deriv_variable

    def _convert_free_variable_deriv(self, eqn, new_variable, cf):
        """
        Create relevant variables/equations when converting a free variable within a derivative.
        :param eqn: the derivative equation containing free variable
        :param new_variable: the new variable representing the converted symbol [new_units]
        :param cf: conversion factor for unit conversion [new units/old units]
        """
        derivative_variable = eqn.args[0].args[0]  # units [x]
        # 1. create a new variable/equation for original derivative
        # will have units [x/old units]
        new_deriv_variable = self._create_new_deriv_variable_and_equation(eqn, derivative_variable)

        # 2. create equation for derivative wrt new variable
        # dx/dnewvar [x/new units] = new_deriv_var [x/old units] / cf [new units/old units]
        expression = sympy.Eq(sympy.Derivative(derivative_variable, new_variable), new_deriv_variable / cf)
        self.add_equation(expression)
        return {'variable': new_deriv_variable, 'expression': eqn.args[0]}

    def _convert_state_variable_deriv(self, original_variable, new_variable, cf):
        """
        Create relevant variables/equations when converting a state variable.
        :param original_variable: the variable to be converted [old units]
        :param new_variable: the new variable representing the converted symbol [new units]
        :param cf: conversion factor for unit conversion [new units/old units]
        :return: a dictionary containing the 'variable' and 'expression' for new derivative
        """
        # 1. find the derivative equation for this variable
        eqn = self._ode_definition_map[original_variable]

        # get free variable symbol
        # units [x]
        wrt_variable = eqn.args[0].args[1]

        # 1. create a new variable/equation for original derivative
        # will have units [old units/x]
        new_deriv_variable = self._create_new_deriv_variable_and_equation(eqn, original_variable)

        # 2. add a new derivative equation
        # dnewvar/dx [new units/x] = new_deriv_var [old units/x] * cf [new units/old units]
        expression = sympy.Eq(sympy.Derivative(new_variable, wrt_variable), new_deriv_variable * cf)
        self.add_equation(expression)
        return {'variable': new_deriv_variable, 'expression': eqn.args[0]}

    def _convert_variable_instance(self, original_variable, cf, units, direction):
        """
        Internal function to create new variable and an equation for it.
        :param original_variable: VariableDummy object to be converted [old units]
        :param cf: conversion factor [new units/old units]
        :param units: Unit object for new units
        :param direction: enumeration value specifying input or output
        :return: the new variable created [new units]
        """
        # Get unique name for new variable
        new_name = self.get_unique_name(original_variable.name + '_converted')

        # If original has initial_value calculate new initial value (only needed for INPUT case)
        new_value = None
        if direction == DataDirectionFlow.INPUT and original_variable.initial_value is not None:
            new_value = original_variable.initial_value * cf

        # Create new variable
        new_variable = self.add_variable(name=new_name, units=units, initial_value=new_value)

        # Transfer cmeta id from original to new variable if the original variable has one
        if original_variable._cmeta_id is not None:
            self.transfer_cmeta_id(original_variable, new_variable)

        # Add/remove/replace equations
        if direction == DataDirectionFlow.INPUT:
            # if direction is input; original var will be replaced by equation so do not need to store initial value
            original_variable.initial_value = None

            # find the equation for the original variable (if any): orig_var = rhs
            # remove equation from model
            # add eqn for new variable in terms of rhs of equation
            #     new_var [new units] = rhs [old units] * cf [new units/old units]
            original_equation = self._var_definition_map.get(original_variable)
            if original_equation is not None:
                new_equation = sympy.Eq(new_variable, original_equation.args[1] * cf)
                self.remove_equation(original_equation)
                self.add_equation(new_equation)

            # add eqn for original variable in terms of new variable
            #     orig_var [old units] = new var [new units] / cf [new units/old units]
            expression = sympy.Eq(original_variable, new_variable / cf)
            self.add_equation(expression, check_duplicates=False)
        else:
            # if direction is output add eqn for new variable in terms of original variable
            #     new_var [new units] = orig_var [old units] * cf [new units/old units]
            expression = sympy.Eq(new_variable, original_variable * cf)
            self.add_equation(expression)

        return new_variable

    def get_unique_name(self, name):
        """
        Creates and returns a unique name, not used in the model.

        :param str name: Suggested unique name.
        :return str: Guaranteed unique name.
        """
        if name in self._name_to_variable:
            name = self.get_unique_name(name + '_a')
        return name


class NumberDummy(sympy.Dummy):
    """
    Used to represent a number with a unit, inside a Sympy expression.

    Unlike sympy expressions, this number type will never be removed in simplify operations etc.

    Number dummies should never be created directly, but always via :meth:`Model.add_number()`

    Assumes the value is real.
    """

    # Sympy annoyingly overwrites __new__
    def __new__(cls, value, *args, **kwargs):
        return super().__new__(cls, str(value), real=True)

    def __init__(self, value, units):
        self.value = float(value)
        self.units = units

    def __float__(self):
        return self.value

    def _eval_evalf(self, prec):
        """This is needed to allow Sympy's ``evalf`` method to represent this value as a float."""
        return sympy.Float(self.value, prec)

    def __str__(self):
        return str(self.value)


class VariableDummy(sympy.Dummy):
    """
    Used to represent a variable (with meta data) in a Sympy expression.

    Variable dummies should never be created directly, but always via :meth:`Model.add_variable()`.

    For the constructor arguments, see :meth:`Model.add_variable()`.

    Assumes the value is real.
    """

    # Sympy annoyingly overwrites __new__
    def __new__(cls, name, *args, **kwargs):
        return super().__new__(cls, name, real=True)

    def __init__(self,
                 name,
                 units,
                 model=None,
                 initial_value=None,
                 public_interface=None,
                 private_interface=None,
                 order_added=None,
                 cmeta_id=None):
        self._model = model
        self.name = name
        self.units = units
        self.initial_value = None if initial_value is None else float(initial_value)

        # Interface properties, only used during parsing
        self.public_interface = public_interface
        self.private_interface = private_interface

        # Variables are either 'source' variables, or receive their value from
        # a variable that they're connected to (using CellML connections).
        # The ``assigned_to`` property is used to indicate where this object
        # receives its value.
        self.assigned_to = None
        if not (private_interface == 'in' or public_interface == 'in'):
            self.assigned_to = self

        # Optional order added, used for sorting sometimes.
        self.order_added = order_added

        # Optional cmeta id
        self._cmeta_id = self._rdf_identity = None
        self._set_cmeta_id(cmeta_id)

        # This variable's type
        # TODO: Define allowed types via enum
        self.type = None

    def __str__(self):
        return self.name

    @property
    def model(self):
        """The :class:`Model` this variable is part of."""
        return self._model

    @property
    def rdf_identity(self):
        """The RDF identity for this variable (will be ``None`` unless the variable has a cmeta id)."""
        return self._rdf_identity

    @property
    def cmeta_id(self):
        """Provides read-only access to the cmeta id."""
        return self._cmeta_id

    def _set_cmeta_id(self, cmeta_id):
        """Sets this variable's cmeta id. Should only be called by Model, which can verify cmeta id uniqueness."""
        self._cmeta_id = cmeta_id
        if cmeta_id is None:
            self._rdf_identity = None
        else:
            self._rdf_identity = create_rdf_node('#' + self._cmeta_id)
