"""This module contains the CellML parser and related classes. It reads a CellML model and stores
model information in the cellmlmanip.Model class. MathML equations are translated to Sympy. RDF is
handled by RDFLib.
"""
import itertools
import os
from collections import OrderedDict, deque
from enum import Enum

from lxml import etree

from cellmlmanip.model import SYMPY_SYMBOL_DELIMITER, Model
from cellmlmanip.transpiler import Transpiler


_UNIT_PREFIXES = {
    'yocto': 1e-24,
    'zepto': 1e-21,
    'atto': 1e-18,
    'femto': 1e-15,
    'pico': 1e-12,
    'nano': 1e-9,
    'micro': 1e-6,
    'milli': 1e-3,
    'centi': 1e-2,
    'deci': 1e-1,
    'deca': 1e+1,
    'deka': 1e+1,
    'hecto': 1e2,
    'kilo': 1e3,
    'mega': 1e6,
    'giga': 1e9,
    'tera': 1e12,
    'peta': 1e15,
    'exa': 1e18,
    'zetta': 1e21,
    'yotta': 1e24
}


class XmlNs(Enum):
    """Namespaces in CellML documents"""
    CELLML = 'http://www.cellml.org/cellml/1.0#'
    CMETA = 'http://www.cellml.org/metadata/1.0#'
    MATHML = 'http://www.w3.org/1998/Math/MathML'
    RDF = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'


class _Component:
    """This hold information about a CellML component. It's for internal-use only. Once the parser
    has created the flattened cellmlmanip.Model instance, components are no longer used"""
    def __init__(self, name):
        self.name = name
        self.parent = None
        self.siblings = set()
        self.encapsulated = set()

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


class Parser(object):
    """Handles parsing of CellML files"""

    @staticmethod
    def with_ns(ns_enum, name):
        """Returns an ElementTree-friendly name with namespace in brackets"""
        return '{%s}%s' % (ns_enum.value, name)

    def __init__(self, filepath):
        """Initialise an instance of Parser

        :param filepath: the full filepath to the CellML model file
        """
        self.filepath = filepath

        # A :class:`Model` object or None
        self.model = None

        # A dictionary mapping component names to _Component objects
        self.components = OrderedDict()

    def parse(self, unit_store=None):
        """
        The main method that reads the XML file and extracts the relevant parts of the CellML model
        definition.

        :param unit_store: Optional :class:`cellmlmanip.units.UnitStore` instance; if given the model will share the
            underlying registry so that conversions between model units and those from the provided store work.
        :return: a :class:`Model` holding CellML model definition, reading for manipulation.
        """

        # Create lxml parser
        parser = etree.XMLParser(no_network=True)

        # Parse, get ElementTree
        tree = etree.parse(self.filepath, parser)

        # Validate CellML syntax
        self._validate(parser, tree)

        # <model> root node - initialise the model object
        model_xml = tree.getroot()
        self.model = Model(model_xml.get('name'),
                           model_xml.get(Parser.with_ns(XmlNs.CMETA, 'id')),
                           unit_store=unit_store)

        # handle the child elements of <model>
        self._add_units(model_xml)
        self._add_rdf(model_xml)

        self._add_components(model_xml)
        self._add_relationships(model_xml)
        self._add_connection(model_xml)

        # Canonicalise representation
        self.model.transform_constants()

        return self.model

    @staticmethod
    def _get_variable_name(component_name, variable_name):
        return component_name + SYMPY_SYMBOL_DELIMITER + variable_name

    def _add_rdf(self, element):
        """
        Finds all ``<RDF>`` definitions under ``<element>`` and adds them to the model.

        :param element: the CellML parent element to search for children RDF tags
        """
        for rdf in element.iter(Parser.with_ns(XmlNs.RDF, 'RDF')):
            self.model.add_rdf(etree.tostring(rdf, encoding=str))

    def _add_units(self, model):
        """
        <model> <units> <unit /> </units> </model>
        :param model: an etree.Element
        """
        units_elements = model.findall(Parser.with_ns(XmlNs.CELLML, 'units'))

        # get list of built-in cellml units
        from cellmlmanip.units import CELLML_UNITS
        units_found = set(CELLML_UNITS)

        # get all the units defined in the cellml model
        definitions_to_add = OrderedDict()
        for units_element in units_elements:
            units_name = units_element.get('name')
            # if it's a defined base unit, we can be add immediately to the model
            if units_element.get('base_units'):
                self.model.units.add_base_unit(units_name)
                units_found.add(units_name)
            # all other units are collected (because they may depend on further user-defined units)
            else:
                unit_elements = [dict(t.attrib) for t in units_element.getchildren()]
                definitions_to_add[units_name] = unit_elements

        iteration = 0
        # while we still have units to add
        while definitions_to_add:
            # get a definition from the top of the list
            unit_name, unit_elements = definitions_to_add.popitem()

            # check whether this unit is defined in terms of units that we know about
            add_now = True
            for unit in unit_elements:
                # if defined in terms of units we don't know about
                if unit['units'] not in units_found:
                    # defer adding this units - add it back to the end of the list
                    definitions_to_add[unit_name] = unit_elements
                    definitions_to_add.move_to_end(unit_name, last=False)
                    add_now = False
                    break

            # unit is defined in terms of known units - ok to add to model
            if add_now:
                definition = self._make_pint_unit_definition(unit_name, unit_elements)
                self.model.units.add_unit(unit_name, definition)
                units_found.add(unit_name)
                iteration = 0
            else:
                # we did not add any units in this iteration - make note
                iteration += 1

                # exit if we have not been able to add any units in the entire list of definitions
                if iteration > len(definitions_to_add):
                    raise ValueError('Cannot create units %s. Cycles or unknown units.' % definitions_to_add)

    def _make_pint_unit_definition(self, units_name, unit_attributes):
        """
        Construct and return a ``pint.UnitDefinition``.

        :param units_name: The unit name
        :param unit_attributes: A list of dicts, where each dict contains the fields
            ``{'multiplier': float, 'units': string, 'exponent': integer, 'prefix': string/integer}``.
            Not all fields are necessary but ``units`` must match a unit in Pint registry.
        """
        full_unit_expr = []

        # For each of the <unit> elements for this unit definition
        for unit_element in unit_attributes:
            # Source the <unit units="XXX"> from our store
            matched_unit = self.model.units.get_unit(unit_element['units'])

            # Construct a string representing the expression for this <unit>
            expr = str(matched_unit)

            # See https://www.cellml.org/specifications/cellml_1.1/#sec_units 5.2.2
            # offset, prefix, exponent, and multiplier

            if 'prefix' in unit_element:
                try:
                    power = _UNIT_PREFIXES[unit_element['prefix']]
                except KeyError:
                    # Assume that prefix is an integer.
                    power = '1e%s' % unit_element['prefix']
                expr = '(%s * %s)' % (expr, power)

            if 'exponent' in unit_element:
                expr = '((%s)**%s)' % (expr, unit_element['exponent'])

            if 'multiplier' in unit_element:
                expr = '(%s * %s)' % (unit_element['multiplier'], expr)

            # Collect/add this particular <unit> definition
            full_unit_expr.append(expr)

        # Join together all the parts of the unit expression and return
        return '*'.join(full_unit_expr)

    def _add_components(self, model):
        """
        <model> <component> </model>
        :param model: an etree.Element
        """
        component_elements = model.findall(Parser.with_ns(XmlNs.CELLML, 'component'))

        # for each component defined in the model
        for element in component_elements:
            # component are only kept in parser to resolve relationships and connections
            name = element.get('name')
            self.components[name] = _Component(name)

            # process the <variable> tags in this component
            variable_to_symbol = self._add_variables(element)

            # process the <math> tags in this component
            self._add_maths(element, variable_to_symbol)

            # Raise error if component units are defined
            component_units = element.findall(Parser.with_ns(XmlNs.CELLML, 'units'))
            if component_units:
                raise ValueError(
                    'Defining units inside components is not supported (found in component ' + name + ').')

            # Raise error if reactions are defined
            reactions = element.findall(Parser.with_ns(XmlNs.CELLML, 'reaction'))
            if reactions:
                raise ValueError(
                    'Reactions are not supported (found in component ' + name + ').')

    def _add_variables(self, component_element):
        """
        <model> <component> <variable> </component> </model>
        :param component_element: an etree.Element
        """
        variable_elements = component_element.findall(Parser.with_ns(XmlNs.CELLML, 'variable'))

        # we keep a {variable name: sympy symbol} lookup that we pass to the transpiler
        variable_lookup_symbol = dict()

        for variable_element in variable_elements:
            attributes = dict(variable_element.attrib)

            # Rename key for cmeta_id (remove namespace from attribute)
            cmeta_id_attribute = Parser.with_ns(XmlNs.CMETA, 'id')
            if cmeta_id_attribute in attributes:
                attributes['cmeta_id'] = attributes.pop(cmeta_id_attribute)

            # mangle the name by prefixing with the component name
            attributes['name'] = Parser._get_variable_name(component_element.get('name'),
                                                           attributes['name'])

            # look up units
            attributes['units'] = self.model.units.get_unit(attributes['units'])

            # model.add_variable() returns sympy dummy created for this variable - keep it
            variable_lookup_symbol[attributes['name']] = self.model.add_variable(**attributes)

        return variable_lookup_symbol

    def _add_maths(self, component_element, variable_to_symbol):
        """
        <model> <component> <math> </component> </model>

        :param component_element: an etree.Element
        :param variable_to_symbol: a ``Dict[str, sympy.Dummy]``
        """
        # get all <math> elements in the component
        math_elements = component_element.findall(Parser.with_ns(XmlNs.MATHML, 'math'))

        # nothing to do if we don't have any <math> elements
        if not math_elements:
            return

        # Method to create symbols
        prefix = component_element.get('name') + SYMPY_SYMBOL_DELIMITER

        def symbol_generator(identifer):
            out = variable_to_symbol.get(prefix + identifer, None)
            assert out is not None, '%s not found in symbol dict' % (prefix + identifer)
            return out

        # reuse transpiler so dummy symbols are kept across <math> elements
        transpiler = Transpiler(
            symbol_generator=symbol_generator,
            number_generator=lambda x, y: self.model.add_number(x, self.model.units.get_unit(y)),
        )

        # for each math element
        for math_element in math_elements:
            sympy_exprs = transpiler.parse_string(etree.tostring(math_element, encoding=str))

            # add each equation from <math> to the model
            for expr in sympy_exprs:
                self.model.add_equation(expr)

    def _add_relationships(self, model: etree.Element):
        group_elements = model.findall(Parser.with_ns(XmlNs.CELLML, 'group'))

        # find all the <group> elements
        for group_element in group_elements:

            # find the relationship for this <group>
            relationship_ref = group_element.findall(Parser.with_ns(XmlNs.CELLML,
                                                                    'relationship_ref'))
            assert len(relationship_ref) == 1
            relationship = relationship_ref[0].attrib.get('relationship')

            # we only handle 'encapsulation' relationships (i.e. ignoring 'containment')
            if relationship == 'encapsulation':
                self._handle_component_ref(group_element, None)

    def _handle_component_ref(self, parent_tag, parent_component):
        # we're going to process all the siblings at the end
        siblings = []

        # for each of the child <component_ref> elements in the parent tag
        for component_ref_element in parent_tag.findall(Parser.with_ns(XmlNs.CELLML,
                                                                       'component_ref')):

            # get the name of the child component
            child_component = component_ref_element.attrib.get('component')

            # add it to the sibling list
            siblings.append(child_component)

            # if we have a parent component for this child component (i.e. not top-level anonymous)
            if parent_component:
                # add the relationship in the component
                self.components[parent_component].add_encapsulated(child_component)
                self.components[child_component].set_parent(parent_component)

            # descend into this <component_ref> tag to handle any children
            self._handle_component_ref(component_ref_element, child_component)

        # if there are siblings in this non-anonymous group
        if parent_component and len(siblings) > 1:
            # register each of the siblings with each other
            for component_a, component_b in itertools.product(siblings, siblings):
                if component_a != component_b:
                    self.components[component_a].add_sibling(component_b)

    def _add_connection(self, model):
        """
        :param model: an etree.Element
        """
        connection_elements = model.findall(Parser.with_ns(XmlNs.CELLML, 'connection'))

        # a list to collect the (source, target) connection tuples
        connect_from_to = []

        # for each connection in the model
        for connection in connection_elements:
            # first child is <map_component>
            child_0 = connection[0]
            assert child_0.tag == Parser.with_ns(XmlNs.CELLML, 'map_components')
            comp_1, comp_2 = (child_0.attrib.get('component_1'),
                              child_0.attrib.get('component_2'))

            # the remaining children are <map_variables> tags
            for child in connection[1:]:
                assert child.tag == Parser.with_ns(XmlNs.CELLML, 'map_variables')
                connect_from_to.append(
                    self._determine_connection_direction(comp_1, child.attrib.get('variable_1'),
                                                         comp_2, child.attrib.get('variable_2'))
                )

        # we add the connection to the model by first connecting
        # those variables we know are source variables
        connections_to_process = deque(connect_from_to)

        # keep processing the list of connections until we've done them all
        unchanged_loop_count = 0
        while connections_to_process:
            # get a connection
            connection = connections_to_process.popleft()

            # if we did not successfully connect this variable (because, e.g., we don't know the
            # source of *this* source)
            if not self.model.connect_variables(*connection):
                # add it back to the list
                connections_to_process.append(connection)
                unchanged_loop_count += 1
            else:
                unchanged_loop_count = 0

            if unchanged_loop_count > len(connections_to_process):
                raise ValueError('Unable to add connections to the model')

    def _determine_connection_direction(self, comp_1, var_1, comp_2, var_2):
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
        """
        def _are_siblings(comp_a, comp_b):
            return self.components[comp_a].parent == self.components[comp_b].parent

        def _parent_of(parent_name, child_name):
            return parent_name == self.components[child_name].parent

        # get the variable information from the model about each end of the connection
        variable_1 = self.model.get_symbol_by_name(self._get_variable_name(comp_1, var_1))
        variable_2 = self.model.get_symbol_by_name(self._get_variable_name(comp_2, var_2))

        # if the components are siblings (either same parent or top-level)
        if _are_siblings(comp_1, comp_2):
            # they are both connected on their public_interface
            if variable_1.public_interface == 'out' and variable_2.public_interface == 'in':
                return variable_1.name, variable_2.name
            elif variable_1.public_interface == 'in' and variable_2.public_interface == 'out':
                return variable_2.name, variable_1.name
        else:
            # determine which component is parent of the other
            if _parent_of(comp_1, comp_2):
                parent_var, child_var = variable_1, variable_2
            else:
                parent_var, child_var = variable_2, variable_1

            # parent/child components are connected using private/public interface, respectively
            if child_var.public_interface == 'in' and parent_var.private_interface == 'out':
                return parent_var.name, child_var.name
            elif child_var.public_interface == 'out' and parent_var.private_interface == 'in':
                return child_var.name, parent_var.name

        raise ValueError('Cannot determine the source & target for connection (%s, %s) - (%s, %s)' %
                         (comp_1, var_1, comp_2, var_2))

    def _validate(self, parser, tree):
        """
        Validates the given lxml ``tree`` against the CellML 1.0 RELAX NG schema.

        :param parser: An `lxml.etree.XMLParser`
        :param tree: An `lxml.etree.ElementTree` made with the given parser.
        """

        # Create RelaxNG object
        path = os.path.join(os.path.dirname(__file__), 'data', 'cellml_1_0.rng')
        rnc = etree.RelaxNG(etree.parse(path, parser))

        # Validate
        if not rnc.validate(tree):
            msg = '. '.join([str(x) for x in rnc.error_log])
            raise ValueError('Invalid or unsupported CellML file. ' + msg)

