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
        self.cellml_model = None

    def parse(self):
        """
        The main method that reads the XML file and extract the relevant parts of CellML model
        definition. Parser class should have been instantiated with the filepath.

        :return: a Model class holding CellML model definition, reading for manipulation
        """
        tree = etree.parse(self.filepath)

        # CellML <model> root node
        model = tree.getroot()
        self.cellml_model = Model(model.get(self.with_ns(XmlNs.CMETA, u'id')))
        for rdf in model.findall(self.with_ns(XmlNs.RDF, u'RDF')):
            self.cellml_model.add_rdf(etree.tostring(rdf, encoding=str))

        # model / units
        self.__add_units(model)

        # model / components
        self.__add_components(model)

        return self.cellml_model

    def __add_units(self, model):
        """  <model> <units> <unit /> </units> </model> """
        units_elements = model.findall(self.with_ns(XmlNs.CML, u'units'))
        for units_element in units_elements:
            units_name = units_element.get(u'name')
            unit_elements = [dict(t.attrib) for t in units_element.getchildren()]
            self.cellml_model.add_unit(units_name, unit_elements)

    def __add_components(self, model):
        """ <model> <component> </model> """
        component_elements = model.findall(self.with_ns(XmlNs.CML, u'component'))
        for component_element in component_elements:
            component_name = component_element.get(u'name')
            self.cellml_model.add_component(component_name)

            # Add all <variable> in this component
            self.__add_variable(component_element)

            # Add any maths for this component
            self.__add_math(component_element)

            # add any RDF for this component
            for rdf in component_element.findall(self.with_ns(XmlNs.RDF, u'RDF')):
                self.cellml_model.add_rdf(etree.tostring(rdf, encoding=str))

    def __add_math(self, component):
        """ <model> <component> <math> </component> </model> """
        component_name = component.get(u'name')

        # NOTE: Only looking for one <math> element
        math_element = component.find(self.with_ns(XmlNs.MATH, u'math'))

        if math_element:
            transpiler = mathml2sympy.Transpiler(dummify=False)
            sympy_exprs = transpiler.parse_string(etree.tostring(math_element, encoding=str))
            self.cellml_model.add_equations(sympy_exprs, component_name)

    def __add_variable(self, component):
        """ <model> <component> <variable> </component> </model> """
        component_name = component.get(u'name')

        variable_elements = component.findall(self.with_ns(XmlNs.CML, u'variable'))
        for variable_element in variable_elements:
            attributes = dict(variable_element.attrib)

            # Rename key for cmeta_id (remove namespace from attribute)
            cmeta_id_attribute = self.with_ns(XmlNs.CMETA, 'id')
            if cmeta_id_attribute in attributes:
                attributes['cmeta_id'] = attributes.pop(cmeta_id_attribute)

            attributes['component_name'] = component_name
            self.cellml_model.add_variable(**attributes)

            # Add any RDF for this <variable>
            for rdf in variable_element.findall(self.with_ns(XmlNs.RDF, u'RDF')):
                self.cellml_model.add_rdf(etree.tostring(rdf, encoding=str))


class Model(object):
    """
    Holds all information about a CellML model and exposes it for manipulation (intention!)
    """
    def __init__(self, identifier):
        self.data = {'cmeta_id': str(identifier), 'units': {}, 'components': [], 'variables': []}
        self.equations = []
        self.rdf = rdflib.Graph()

    def add_unit(self, units_name: str, unit_elements: dict):
        """
        Adds information about <units> in <model>
        """
        self.data['units'][units_name] = unit_elements

    def add_component(self, name: str):
        """
        Adds name to list of <component>s in the <model>
        """
        self.data['components'].append(name)

    def add_variable(self, name, component_name, units,
                     initial_value=None, private_interface=None,
                     public_interface=None, cmeta_id=None):
        """
        Adds a <variable> of <component> of <model>
        """
        # Create a dict from the method arguments but don't store None-s or self
        # TODO: better way to do this?
        attributes = {key: value for (key, value) in locals().items()
                      if value is not None and key is not 'self'}
        self.data['variables'].append(attributes)

    def add_rdf(self, rdf: str):
        """
        Takes RDF string and stores it in an RDFlib.Graph for the model. Can be called repeatedly.
        """
        self.rdf.parse(StringIO(rdf))

    def add_equations(self, sympy_exprs, component_name):
        """
        Adds <math> equation(s) for <component> for <model>. The equations are Sympy equality
        expressions
        """
        for expr in sympy_exprs:
            self.equations.append({'component_name': component_name, 'eqn': expr})

    def __str__(self):
        # Can be removed later - here for development
        class _CustomEncoder(json.JSONEncoder):
            """Converts Sympy expressions and RDF graphs to str"""
            def default(self, o):
                if isinstance(o, sympy.Eq):
                    return str(o)
                elif isinstance(o, rdflib.Graph):
                    return [str(t) for t in o]
                return json.JSONEncoder.default(self, o)

        everything = {'data': self.data, 'equations': self.equations, 'rdf': self.rdf}
        return json.dumps(everything, indent=4, cls=_CustomEncoder)
