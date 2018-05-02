"""
This module contains the CellML parser. It reads CellML model and stores model information in a
CellML Model class. MathML equations are translated to Sympy. RDF is handled by RDFLib.
"""
from enum import Enum
import json
from io import StringIO
import logging
import sys

from lxml import etree
import rdflib
import sympy

from cellmlmanip import mathml2sympy

logging.basicConfig(level=logging.DEBUG,
                    format="%(name)s: %(levelname)s: %(message)s",
                    stream=sys.stderr)
logging.getLogger().handlers[0].setLevel(logging.DEBUG)


class XmlNs(Enum):
    """
    Standard namespaces present in CellML documents
    TODO: Other common namespaces?
    TODO: Should we be picking these up from somewhere else?
    """
    CML = 'http://www.cellml.org/cellml/1.0#'
    CMETA = 'http://www.cellml.org/metadata/1.0#'
    MATH = 'http://www.w3.org/1998/Math/MathML'
    RDF = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'


class Parser(object):
    """
    This class handles parsing of CellML files. It is instantiated with the CellML filepath.
    """
    @staticmethod
    def with_ns(ns_enum, tag):
        """Returns an ElementTree friendly tag with namespace in brackets"""
        return u'{%s}%s' % (ns_enum.value, tag)

    def __init__(self, filepath):
        self.filepath = filepath
        self.model = None

    def parse(self):
        """
        The main method that reads the XML file and extract the relevant parts of CellML model
        definition. Parser class should have been instantiated with the filepath.

        :return: a Model class holding CellML model definition, reading for manipulation
        """
        tree = etree.parse(self.filepath)

        # <model> root node - initialise the model object
        model_xml = tree.getroot()
        self.model = Model(model_xml.get(self.with_ns(XmlNs.CMETA, u'id')))

        # handle the child elements of <model>
        self.__add_units(model_xml)
        self.__add_components(model_xml)
        self.__add_rdf(model_xml)

        return self.model

    def __add_rdf(self, element):
        for rdf in element.findall(self.with_ns(XmlNs.RDF, u'RDF')):
            self.model.add_rdf(etree.tostring(rdf, encoding=str))

    def __add_units(self, model):
        """  <model> <units> <unit /> </units> </model> """
        units_elements = model.findall(self.with_ns(XmlNs.CML, u'units'))
        for units_element in units_elements:
            units_name = units_element.get(u'name')
            unit_elements = [dict(t.attrib) for t in units_element.getchildren()]
            self.model.add_unit(units_name, unit_elements)

    def __add_components(self, model):
        """ <model> <component> </model> """
        component_elements = model.findall(self.with_ns(XmlNs.CML, u'component'))

        # for each component defined in the model
        for component_element in component_elements:
            # create an instance of Component
            component = Component(component_element.get(u'name'))

            # Add the child elements under <component>
            self.__add_variables(component, component_element)
            self.__add_math(component, component_element)
            self.__add_rdf(component_element)

            # Add the component instance to the model
            self.model.add_component(component)

    def __add_math(self, component, component_element):
        """ <model> <component> <math> </component> </model> """
        # NOTE: Only looking for one <math> element
        math_element = component_element.find(self.with_ns(XmlNs.MATH, u'math'))
        if math_element is not None:
            transpiler = mathml2sympy.Transpiler(dummify=True)
            # TODO: check whether element can be passed directly without .tostring()
            sympy_exprs = transpiler.parse_string(etree.tostring(math_element, encoding=str))
            component.equations = sympy_exprs
            component.add_symbol_info(transpiler.metadata)

    def __add_variables(self, component, component_element):
        """ <model> <component> <variable> </component> </model> """
        variable_elements = component_element.findall(self.with_ns(XmlNs.CML, u'variable'))
        for variable_element in variable_elements:
            attributes = dict(variable_element.attrib)

            # Rename key for cmeta_id (remove namespace from attribute)
            cmeta_id_attribute = self.with_ns(XmlNs.CMETA, 'id')
            if cmeta_id_attribute in attributes:
                attributes['cmeta:id'] = attributes.pop(cmeta_id_attribute)

            component.variables[attributes['name']] = attributes

            # Add any RDF for this <variable>
            self.__add_rdf(variable_element)


class Component(object):

    def __init__(self, name):
        self.name = name
        self.variables = dict()
        self.equations = None
        self.numbers = {}

    def __str__(self):
        return "%s\n\tequations: %s\n\tvariables: %s\n\tnumbers: %s" % (self.name, self.equations, self.variables, self.numbers)

    def add_symbol_info(self, metadata):
        """
        Called after adding equations and variables from the CellML model
        Collects the disparate bits of information about the equation into a single place
        """
        if not self.equations:
            return

        # Get all unique symbols in the equations
        dummy_symbols = set().union(*[e.free_symbols for e in self.equations])

        # For each symbol in these equations
        for s in dummy_symbols:
            # if the symbol is one that's defined as a <variable> in the component
            if s.name in self.variables:
                self.variables[s.name]['sympy.Dummy'] = s
                if s in metadata:
                    self.variables[s.name].update(metadata[s])

            # create variable entry for dummified derivative
            if s in metadata:
                if 'sympy.Derivative' in metadata[s]:
                    derivative = metadata[s]['sympy.Derivative']
                    y_symbol = derivative.free_symbols.pop()
                    x_symbol = derivative.variables[0]
                    # the bound and wrt symbols should be <variable>s in the component
                    self.variables[y_symbol.name]['sympy.Dummy'] = y_symbol
                    self.variables[x_symbol.name]['sympy.Dummy'] = x_symbol
                    self.variables[s.name] = metadata[s]
                elif 'sympy.Number' in metadata[s]:
                    self.numbers[s] = metadata[s]


class Model(object):
    """
    Holds all information about a CellML model and exposes it for manipulation (intention!)
    """
    def __init__(self, identifier):
        self.name = identifier
        self.units = {}
        self.components = {}
        self.rdf = rdflib.Graph()

    def add_unit(self, units_name: str, unit_elements: dict):
        """
        Adds information about <units> in <model>
        """
        self.units[units_name] = unit_elements

    def add_component(self, component: Component):
        """
        Adds name to list of <component>s in the <model>
        """
        self.components[component.name] = component

    def add_rdf(self, rdf: str):
        """
        Takes RDF string and stores it in an RDFlib.Graph for the model. Can be called repeatedly.
        """
        self.rdf.parse(StringIO(rdf))
