"""
This module contains the CellML parser. It reads CellML model and stores model information in a
CellML Model class. MathML equations are translated to Sympy. RDF is handled by RDFLib.
"""
from enum import Enum
from typing import Dict

from lxml import etree

from cellmlmanip import mathml2sympy
from cellmlmanip.model import Component, Model


class XmlNs(Enum):
    """
    Standard namespaces present in CellML documents
    TODO: Other common namespaces?
    TODO: Should we be picking these up from somewhere else?
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
        return u'{%s}%s' % (ns_enum.value, name)

    def __init__(self, filepath: str) -> None:
        """Initialise an instance of Parser

        :param filepath: the full filepath to the CellML model file
        """
        self.filepath: str = filepath
        self.model: Model = None

    def parse(self):
        """The main method that reads the XML file and extract the relevant parts of CellML model
        definition. Parser class should have been instantiated with the filepath.

        :return: a Model class holding CellML model definition, reading for manipulation
        """
        tree = etree.parse(self.filepath)

        # <model> root node - initialise the model object
        model_xml = tree.getroot()
        self.model = Model(model_xml.get(Parser.with_ns(XmlNs.CMETA, u'id')))

        # handle the child elements of <model>
        self.__add_units(model_xml)
        self.__add_components(model_xml)
        self.__add_rdf(model_xml)
        self.__add_connection(model_xml)

        return self.model

    def __add_rdf(self, element: etree.Element):
        """Finds all <RDF> definitions under <element> and adds them to the model

        :param element: the CellML parent element to search for children RDF tags
        """
        for rdf in element.findall(Parser.with_ns(XmlNs.RDF, u'RDF')):
            self.model.add_rdf(etree.tostring(rdf, encoding=str))

    def __add_units(self, model: etree.Element):
        """  <model> <units> <unit /> </units> </model> """
        units_elements = model.findall(Parser.with_ns(XmlNs.CELLML, u'units'))
        units_collected = {}
        for units_element in units_elements:
            units_name = units_element.get(u'name')
            unit_elements = [dict(t.attrib) for t in units_element.getchildren()]
            units_collected[units_name] = unit_elements
        self.model.add_unit(units_collected)

    def __add_components(self, model: etree.Element):
        """ <model> <component> </model> """
        component_elements = model.findall(Parser.with_ns(XmlNs.CELLML, u'component'))

        # for each component defined in the model
        for component_element in component_elements:
            # create an instance of Component
            component = Component(component_element.get(u'name'), self.model)

            # Add the child elements under <component>
            self.__add_variables(component, component_element)
            self.__add_maths(component, component_element)
            self.__add_rdf(component_element)

            # Add the component instance to the model
            self.model.add_component(component)

    def __add_maths(self, component: Component, component_element: etree.Element):
        """ <model> <component> <math> </component> </model> """
        # NOTE: Only looking for one <math> element
        math_element = component_element.find(Parser.with_ns(XmlNs.MATHML, u'math'))
        if math_element is not None:
            transpiler = mathml2sympy.Transpiler(dummify=True, symbol_prefix=component.name + '__')
            # TODO: check whether element can be passed directly without .tostring()
            sympy_exprs = transpiler.parse_string(etree.tostring(math_element, encoding=str))
            component.equations.extend(sympy_exprs)
            component.collect_variable_attributes(transpiler.metadata)

    def __add_variables(self, component: Component, component_element: etree.Element):
        """ <model> <component> <variable> </component> </model> """
        variable_elements = component_element.findall(Parser.with_ns(XmlNs.CELLML, u'variable'))
        for variable_element in variable_elements:
            attributes: Dict = dict(variable_element.attrib)

            # Rename key for cmeta_id (remove namespace from attribute)
            cmeta_id_attribute = Parser.with_ns(XmlNs.CMETA, 'id')
            if cmeta_id_attribute in attributes:
                attributes['cmeta:id'] = attributes.pop(cmeta_id_attribute)

            component.variables[attributes['name']] = attributes

            # Add any RDF for this <variable>
            self.__add_rdf(variable_element)

    def __add_connection(self, model: etree.Element):
        connection_elements = model.findall(Parser.with_ns(XmlNs.CELLML, u'connection'))
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
