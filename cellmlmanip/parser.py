"""
This module contains the CellML parser. It reads CellML model and stores model information in a
CellML Model class. MathML equations are translated to Sympy. RDF is handled by RDFLib.
"""
import itertools
from collections import OrderedDict, deque
from enum import Enum
from typing import Dict

from lxml import etree

from cellmlmanip import mathml2sympy
from cellmlmanip.model import SYMPY_SYMBOL_DELIMITER, Model, Variable


class XmlNs(Enum):
    """
    Standard namespaces present in CellML documents
    """
    CELLML = 'http://www.cellml.org/cellml/1.0#'
    CMETA = 'http://www.cellml.org/metadata/1.0#'
    MATHML = 'http://www.w3.org/1998/Math/MathML'
    RDF = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'


class ComponentNew:
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

    def __str__(self) -> str:
        return '%s(%s)' % (
            type(self).__name__,
            ', '.join('%s=%s' % item for item in vars(self).items() if item[1])
        )


class Parser(object):
    """Handles parsing of CellML files
    """

    @staticmethod
    def with_ns(ns_enum, name):
        """Returns an ElementTree-friendly name with namespace in brackets"""
        return '{%s}%s' % (ns_enum.value, name)

    def __init__(self, filepath: str) -> None:
        """Initialise an instance of Parser

        :param filepath: the full filepath to the CellML model file
        """
        self.filepath: str = filepath
        self.model: Model = None
        self.components: Dict[str, ComponentNew] = OrderedDict()

    def parse(self) -> Model:
        """The main method that reads the XML file and extract the relevant parts of CellML model
        definition. Parser class should have been instantiated with the filepath.

        :return: a Model class holding CellML model definition, reading for manipulation
        """
        tree = etree.parse(self.filepath)

        # <model> root node - initialise the model object
        model_xml = tree.getroot()
        self.model = Model(model_xml.get(Parser.with_ns(XmlNs.CMETA, 'id')))

        # handle the child elements of <model>
        self._add_units(model_xml)
        self._add_rdf(model_xml)

        self._add_components(model_xml)
        self._add_relationships(model_xml)
        self._add_connection(model_xml)

        return self.model

    @staticmethod
    def _get_variable_name(component_name, variable_name):
        return component_name + SYMPY_SYMBOL_DELIMITER + variable_name

    def _add_rdf(self, element: etree.Element):
        """Finds all <RDF> definitions under <element> and adds them to the model

        :param element: the CellML parent element to search for children RDF tags
        """
        for rdf in element.iter(Parser.with_ns(XmlNs.RDF, 'RDF')):
            self.model.add_rdf(etree.tostring(rdf, encoding=str))

    def _add_units(self, model: etree.Element):
        """  <model> <units> <unit /> </units> </model> """
        units_elements = model.findall(Parser.with_ns(XmlNs.CELLML, 'units'))
        for units_element in units_elements:
            units_name = units_element.get('name')
            if units_element.get('base_units'):
                self.model.add_unit(units_name, unit_attributes=None, base_units=True)
            else:
                unit_elements = [dict(t.attrib) for t in units_element.getchildren()]
                self.model.add_unit(units_name, unit_attributes=unit_elements)

    def _add_components(self, model: etree.Element):
        """ <model> <component> </model> """
        component_elements = model.findall(Parser.with_ns(XmlNs.CELLML, 'component'))

        # for each component defined in the model
        for component_element in component_elements:
            # create an instance of Component
            self.components[component_element.get('name')] = ComponentNew(component_element.get('name'))

            # Add the child elements under <component>
            variables = self._add_variables(component_element)
            # print('variables:\n', variables)
            self._add_maths(component_element, variables)

    def _add_maths(self, component_element: etree.Element, symbol_lookup):
        """ <model> <component> <math> </component> </model> """
        # get all <math> elements in the component
        math_elements = component_element.findall(Parser.with_ns(XmlNs.MATHML, 'math'))

        # nothing to do if we don't have any <math> elements
        if not math_elements:
            return

        # reuse transpiler so cache of dummy symbols is preserved across <math> elements
        transpiler = mathml2sympy.Transpiler(
            dummify=True,
            symbol_prefix=component_element.get('name') + SYMPY_SYMBOL_DELIMITER,
            symbol_lookup=symbol_lookup
        )

        # for each math element
        for math_element in math_elements:
            # TODO: check whether element can be passed directly without .tostring()
            sympy_exprs = transpiler.parse_string(etree.tostring(math_element, encoding=str))
            for expr in sympy_exprs:
                self.model.add_equation(expr)

        # Add metadata collected whilst parsing <math> elements to the model
        # if transpiler.metadata:
        #     print(transpiler.metadata)
        for symbol, attributes in transpiler.metadata.items():
            self.model.add_number(symbol, attributes)

    def _add_variables(self, component_element: etree.Element):
        """ <model> <component> <variable> </component> </model> """
        variable_elements = component_element.findall(Parser.with_ns(XmlNs.CELLML, 'variable'))
        variable_lookup_symbol = dict()
        for variable_element in variable_elements:
            attributes: Dict = dict(variable_element.attrib)

            # Rename key for cmeta_id (remove namespace from attribute)
            cmeta_id_attribute = Parser.with_ns(XmlNs.CMETA, 'id')
            if cmeta_id_attribute in attributes:
                attributes['cmeta_id'] = attributes.pop(cmeta_id_attribute)

            attributes['_original_name'] = attributes['name']
            attributes['_component_name'] = component_element.get('name')

            # mangle the name by prefixing with the component name
            attributes['name'] = Parser._get_variable_name(component_element.get('name'), attributes['name'])
            variable_lookup_symbol[attributes['name']] = self.model.add_variable(**attributes)
        return variable_lookup_symbol

    def _add_connection(self, model: etree.Element):
        connection_elements = model.findall(Parser.with_ns(XmlNs.CELLML, 'connection'))
        source_target_connections = []
        for connection in connection_elements:
            # there can be one <map_component> and multiple <map_variables>
            map_variables = []
            for child in connection:
                # There should *always* be at least one <map_component>
                if child.tag == Parser.with_ns(XmlNs.CELLML, 'map_components'):
                    map_component = (child.attrib.get('component_1'),
                                     child.attrib.get('component_2'))
                elif child.tag == Parser.with_ns(XmlNs.CELLML, 'map_variables'):
                    map_variables.append((child.attrib.get('variable_1'),
                                          child.attrib.get('variable_2')))
            for variable_0, variable_1 in map_variables:
                source_target_connections.append(self._connect(map_component[0], variable_0, map_component[1], variable_1))

        # confirm assignment of those variables that are not connected to anything else
        self.model.assigned_unconnected_variables()

        # make the connections in the model
        connections_to_process = deque(source_target_connections)

        # For testing: shuffle the order of connections
        # TODO: REMOVE!
        from random import shuffle
        shuffle(connections_to_process)

        # While we still have connections left to process
        while connections_to_process:
            # Get connection at front of queue
            connection = connections_to_process.popleft()
            success = self.model.connect_variables(*connection)
            if not success:
                connections_to_process.append(connection)
            # TODO: track looping and break if we can't exit

    def _connect(self, comp_1, var_1: Variable, comp_2, var_2: Variable):
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

        n1 = self._get_variable_name(comp_1, var_1)
        n2 = self._get_variable_name(comp_2, var_2)

        V1 = self.model.variables_x[n1]
        V2 = self.model.variables_x[n2]

        def _are_siblings(comp_a, comp_b):
            return self.components[comp_a].parent == self.components[comp_b].parent

        def _parent_of(parent_name, child_name):
            return parent_name == self.components[child_name].parent

        # if the components are siblings (either same parent or top-level)
        if _are_siblings(comp_1, comp_2):
            # they are both connected on their public_interface
            if V1.public_interface == 'out' and V2.public_interface == 'in':
                return (n1, n2)
            elif V1.public_interface == 'in' and V2.public_interface == 'out':
                return (n2, n1)
        else:
            # determine which component is parent of the other
            if _parent_of(comp_1, comp_2):
                parent_comp, child_comp = comp_1, comp_2
                parent_var, child_var = var_1, var_2
                parent_V, child_V = V1, V2
            else:
                parent_comp, child_comp = comp_2, comp_1
                parent_var, child_var = var_2, var_1
                parent_V, child_V = V2, V1

            # parent/child components are connected using private/public interface, respectively
            if child_V.public_interface == 'in' and parent_V.private_interface == 'out':
                return (self._get_variable_name(parent_comp, parent_var),
                        self._get_variable_name(child_comp, child_var))
            elif child_V.public_interface == 'out' and parent_V.private_interface == 'in':
                return (self._get_variable_name(child_comp, child_var),
                        self._get_variable_name(parent_comp, parent_var))

        raise ValueError('Cannot determine the source & target for connection (%s, %s) - (%s, %s)' %
                         (comp_1, var_1, comp_2, var_2))

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

