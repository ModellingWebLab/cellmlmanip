import logging
import os

import pytest
import sympy

from cellmlmanip import parser

logging.getLogger().setLevel(logging.INFO)


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
        assert len(model.units) == 11  # grep -c '<units ' test_simple_odes.cellml

    def test_component_count(self, model):
        assert len(model.components) == 17  # grep -c '<component ' test_simple_odes.cellml

    def test_equations_count(self, model):
        equation_count = 0
        for component in model.components.values():
            equation_count += len(component.equations)
        assert equation_count == 16  # determined by hand

    def test_variable_find(self, model):
        assert model.find_variable({'cmeta:id': 'time'}) == [{'cmeta:id': 'time',
                                                              'name': 'time',
                                                              'public_interface': 'out',
                                                              'units': 'ms'}]
        matched = model.find_variable({'cmeta:id': 'sv12'})
        assert len(matched) == 1 and matched[0]['sympy.Dummy'].name == 'sv1'

    def test_connections_loaded(self, model):
        assert len(model.connections) == 26  # grep -c '<map_variables ' test_simple_odes.cellml
        first, second = model.connections[0]
        component_1, variable_1 = first[0], first[1]
        component_2, variable_2 = second[0], second[1]
        var_one = model.components[component_1].variables[variable_1]
        var_two = model.components[component_2].variables[variable_2]
        assert component_1 == 'single_independent_ode'
        assert var_one['name'] == 'time' and var_one['public_interface'] == 'in'
        assert component_2 == 'environment'
        assert var_two['name'] == 'time' and var_two['public_interface'] == 'out'

    def test_connections(self, model):
        model.make_connections()
        # Check environment component's time variable has propagated
        environment__time = model.components['environment'].variables['time']['assignment']

        # We're checking sympy.Dummy objects (same name != same hash)
        assert isinstance(environment__time, sympy.Dummy)
        assert environment__time != sympy.Dummy('time')

        state_units_conversion2__time = \
            model.components['state_units_conversion2'].variables['time']['assignment']
        assert environment__time == state_units_conversion2__time

        # propagated environment time to inside nested component circle_y
        circle_y__time = model.components['circle_y'].variables['time']['assignment']
        assert environment__time == circle_y__time

        # we have a new equation that links together times in different units
        time_units_conversion2__time = \
            model.components['time_units_conversion2'].variables['time']['assignment']
        equation = sympy.Eq(time_units_conversion2__time, environment__time)
        assert equation in model.components['time_units_conversion2'].equations
