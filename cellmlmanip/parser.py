"""
This module contains the CellML parser. It reads CellML model and stores model information in a
CellML Model class. MathML equations are translated to Sympy. RDF is handled by RDFLib.
"""
import itertools
from enum import Enum
from typing import Dict

from lxml import etree

from cellmlmanip import mathml2sympy
from cellmlmanip.model import SYMPY_SYMBOL_DELIMITER, Component, Model


class XmlNs(Enum):
    """
    Standard namespaces present in CellML documents
    """
    CELLML = 'http://www.cellml.org/cellml/1.0#'
    CMETA = 'http://www.cellml.org/metadata/1.0#'
    MATHML = 'http://www.w3.org/1998/Math/MathML'
    RDF = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'


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
        self._add_components(model_xml)
        self._add_relationships(model_xml)
        self._add_rdf(model_xml)
        self._add_connection(model_xml)

        return self.model

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
            component = Component(component_element.get('name'), self.model)

            # Add the child elements under <component>
            self._add_variables(component, component_element)
            self._add_maths(component, component_element)

            # Add the component instance to the model
            self.model.add_component(component)

    def _add_maths(self, component: Component, component_element: etree.Element):
        """ <model> <component> <math> </component> </model> """
        # get all <math> elements in the component
        math_elements = component_element.findall(Parser.with_ns(XmlNs.MATHML, 'math'))

        # nothing to do if we don't have any <math> elements
        if not math_elements:
            return

        # reuse transpiler so cache of dummy symbols is preserved across <math> elements
        transpiler = mathml2sympy.Transpiler(
            dummify=True, symbol_prefix=component.name+SYMPY_SYMBOL_DELIMITER,
        )

        # for each math element
        for math_element in math_elements:
            # TODO: check whether element can be passed directly without .tostring()
            sympy_exprs = transpiler.parse_string(etree.tostring(math_element, encoding=str))
            component.equations.extend(sympy_exprs)

        # Add metadata collected whilst parsing <math> elements to the model
        component.collect_variable_attributes(transpiler.metadata)

    def _add_variables(self, component: Component, component_element: etree.Element):
        """ <model> <component> <variable> </component> </model> """
        variable_elements = component_element.findall(Parser.with_ns(XmlNs.CELLML, 'variable'))
        for variable_element in variable_elements:
            attributes: Dict = dict(variable_element.attrib)

            # Rename key for cmeta_id (remove namespace from attribute)
            cmeta_id_attribute = Parser.with_ns(XmlNs.CMETA, 'id')
            if cmeta_id_attribute in attributes:
                attributes['cmeta:id'] = attributes.pop(cmeta_id_attribute)

            component.variables[attributes['name']] = attributes

    def _add_connection(self, model: etree.Element):
        connection_elements = model.findall(Parser.with_ns(XmlNs.CELLML, 'connection'))
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
                self.model.connections.append(((map_component[0], variable_0),
                                               (map_component[1], variable_1)))

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
                self.model.components[parent_component].add_encapsulated(child_component)
                self.model.components[child_component].set_parent(parent_component)

            # descend into this <component_ref> tag to handle any children
            self._handle_component_ref(component_ref_element, child_component)

        # if there are siblings in this non-anonymous group
        if parent_component and len(siblings) > 1:
            # register each of the siblings with each other
            for component_a, component_b in itertools.product(siblings, siblings):
                if component_a != component_b:
                    self.model.components[component_a].add_sibling(component_b)
