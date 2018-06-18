import logging
import os

import pytest
import sympy
from cellmlmanip.model import QuantityStore
from sympy.physics.units import Quantity

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
        assert len(matched) == 1 and matched[0]['sympy.Dummy'].name == 'single_ode_rhs_const_var__sv1'

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
        assert environment__time != sympy.Dummy(environment__time.name)

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

    def test_quantity_translation(self, model):
        import sympy.physics.units as u

        assert model.units.get_quantity('dimensionless').scale_factor == sympy.Rational(1, 1)

        # Units defined in the test CellML <model>:
        unit_names = ['ms', 'per_ms', 'usec', 'mV', 'per_mV', 'uV', 'mV_per_ms', 'mV_per_s',
                      'mV_per_usec', 'mM', 'mM_per_ms']
        for name in unit_names:
            model.units.get_quantity(name)

        # Sympy built-in units
        assert model.units.get_quantity('ms') == u.millisecond
        assert model.units.get_quantity('centimeter') == u.centimeter

        # Custom units defined in CellML example
        assert u.convert_to(model.units.get_quantity('per_ms'), u.millisecond) == 1/u.millisecond
        assert u.convert_to(model.units.get_quantity('usec'), u.microsecond) == u.microsecond
        assert u.convert_to(
            model.units.get_quantity('mM_per_ms'),
            [u.mole, u.liter, u.millisecond]) == (u.mole / 1000) / (u.liter * u.millisecond)

    def test_add_units_to_equations(self, model):
        for _, component in model.components.items():
            component.add_units_to_equations()

            for equation in component.equations:
                print(equation)
                lhs_units = QuantityStore.summarise_units(equation.lhs)
                rhs_units = QuantityStore.summarise_units(equation.rhs)
                if QuantityStore.is_equal(lhs_units, rhs_units):
                    print('\t', lhs_units, '==', rhs_units)
                else:
                    print('\t', lhs_units, '!=', rhs_units)

        # mV/millisecond == mV_per_ms
        test_equation = model.components['single_independent_ode'].equations[0]
        lhs_units = QuantityStore.summarise_units(test_equation.lhs)
        rhs_units = QuantityStore.summarise_units(test_equation.rhs)
        assert QuantityStore.is_equal(lhs_units, rhs_units)

        # mV_per_usec != uV/millisecond
        test_equation = model.components['deriv_on_rhs2b'].equations[0]
        lhs_units = QuantityStore.summarise_units(test_equation.lhs)
        rhs_units = QuantityStore.summarise_units(test_equation.rhs)
        assert not QuantityStore.is_equal(lhs_units, rhs_units)

    def test_unit_extraction(self):
        import sympy.physics.units as units

        eq = (5*units.mile/(2*units.hour + 10*units.minute))**(8*units.gram)
        assert QuantityStore.summarise_units(eq) == (units.mile/(units.hour + units.minute))**units.gram

        millivolts = Quantity('millivolts', units.voltage, units.milli * units.volts, 'mV')
        x, y = sympy.symbols('x y')
        eq = (millivolts / units.millisecond)*sympy.Derivative(x, y)
        assert QuantityStore.summarise_units(eq) == (millivolts / units.milliseconds)

    def _test_print(self, model):
        # show equations
        for _, component in model.components.items():
            for equation in component.equations:
                print(equation)
