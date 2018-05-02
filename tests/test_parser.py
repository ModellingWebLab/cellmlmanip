import os
import logging

import pytest


from cellmlmanip import parser

logging.basicConfig()
logging.getLogger().setLevel(logging.WARNING)


class TestParser(object):
    @staticmethod
    def get_test_cellml_filepath():
        return os.path.join(os.path.dirname(__file__), "cellml_files", "test_simple_odes.cellml")

    @pytest.fixture(scope="class")
    def model(self):
        p = parser.Parser(TestParser.get_test_cellml_filepath())
        model = p.parse()
        return model

    def test_unit_count(self, model):
        assert len(model.units) == 11

    def test_component_count(self, model):
        assert len(model.components) == 17

    def test_equations_count(self, model):
        equation_count = 0
        for component in model.components.values():
            equation_count += len(component.equations)
        assert equation_count == 16

    def test_variable_find(self, model):
        assert model.find_variable({'cmeta:id': 'time'}) == [{'cmeta:id': 'time',
                                                              'name': 'time',
                                                              'public_interface': 'out',
                                                              'units': 'ms'}]
        matched = model.find_variable({'cmeta:id': 'sv12'})
        assert len(matched) == 1 and matched[0]['sympy.Dummy'].name == 'sv1'

    def test_connections(self, model):
        assert len(model.connections) == 26
        first, second = model.connections[0]
        component_1, variable_1 = first[0], first[1]
        component_2, variable_2 = second[0], second[1]
        var_one = model.components[component_1].variables[variable_1]
        var_two = model.components[component_2].variables[variable_2]
        assert component_1 == 'single_independent_ode'
        assert var_one['name'] == 'time' and var_one['public_interface'] == 'in'
        assert component_2 == 'environment'
        assert var_two['name'] == 'time' and var_two['public_interface'] == 'out'

