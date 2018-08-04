import logging
import os

import pytest
import sympy
import sympy.physics.units as units

from cellmlmanip import parser
from cellmlmanip.model import QuantityStore

logging.getLogger().setLevel(logging.INFO)


class TestParser(object):

    @pytest.fixture(scope="class")
    def model(self):
        """Parses example CellML and returns model"""
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "test_simple_odes.cellml"
        )
        p = parser.Parser(example_cellml)
        model = p.parse()
        return model

    def test_component_count(self, model):
        assert len(model.components) == 21  # grep -c '<component ' test_simple_odes.cellml

    def test_equations_count(self, model):
        equation_count = 0
        for component in model.components.values():
            equation_count += len(component.equations)
        assert equation_count == 18  # determined by hand

    def test_variable_find(self, model):
        assert model.find_variable({'cmeta:id': 'time'}) == [{'cmeta:id': 'time',
                                                              'name': 'time',
                                                              'public_interface': 'out',
                                                              'units': 'ms'}]
        matched = model.find_variable({'cmeta:id': 'sv12'})
        assert len(matched) == 1 and \
            matched[0]['sympy.Dummy'].name == 'single_ode_rhs_const_var__sv1'

    def test_connections_loaded(self, model):
        assert len(model.connections) == 32  # grep -c '<map_variables ' test_simple_odes.cellml
        first, second = model.connections[0]
        component_1, variable_1 = first
        component_2, variable_2 = second
        var_one = model.components[component_1].variables[variable_1]
        var_two = model.components[component_2].variables[variable_2]
        assert component_1 == 'single_independent_ode'
        assert var_one['name'] == 'time'
        assert var_one['public_interface'] == 'in'
        assert component_2 == 'environment'
        assert var_two['name'] == 'time'
        assert var_two['public_interface'] == 'out'

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
        assert model.units.get_quantity('dimensionless').scale_factor == sympy.Rational(1, 1)

        # Units defined in the test CellML <model>:
        unit_names = ['ms', 'per_ms', 'usec', 'mV', 'per_mV', 'uV', 'mV_per_ms', 'mV_per_s',
                      'mV_per_usec', 'mM', 'mM_per_ms']
        for name in unit_names:
            model.units.get_quantity(name)

        # Sympy built-in units
        assert model.units.get_quantity('kilogram') == units.kilogram

        # Custom units defined in CellML example
        assert units.convert_to(model.units.get_quantity('per_ms'), 1/units.millisecond) == 1/units.millisecond
        assert units.convert_to(model.units.get_quantity('usec'), units.microsecond) == units.microsecond
        assert units.convert_to(
            model.units.get_quantity('mM_per_ms'),
            [units.mole, units.liter, units.millisecond]) == (units.mole / 1000) / (units.liter * units.millisecond)

    def test_add_units_to_equations(self, model):
        # This is an irreversible operation # TODO: don't mutate?
        for component in model.components.values():
            component.add_units_to_equations()

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

        # Check a specific RHS->LHS unit conversion
        test_equation = model.components['deriv_on_rhs2b'].equations[0]
        new_rhs = units.convert_to(test_equation.rhs, lhs_units)
        new_rhs_units = QuantityStore.summarise_units(new_rhs)
        assert QuantityStore.is_equal(lhs_units, new_rhs_units)

        # TODO: work in progress...trying to understand what's going on here
        def simplify_units_until_no_change(expr):
            current_expression = expr
            while True:
                new_expression = QuantityStore.summarise_units(current_expression)
                if current_expression == new_expression:
                    break
                current_expression = new_expression
            return new_expression

        # Try fixing all units on the RHS so that they match the LHS
        for component in model.components.values():
            for index, equation in enumerate(component.equations):
                lhs_units = QuantityStore.summarise_units(equation.lhs)
                rhs_units = simplify_units_until_no_change(equation.rhs)
                if not QuantityStore.is_equal(lhs_units, rhs_units):
                    new_rhs = units.convert_to(equation.rhs, lhs_units)
                    # Create a new equality with the converted RHS and replace original
                    equation = sympy.Eq(equation.lhs, new_rhs)
                    component.equations[index] = equation
                    lhs_units = QuantityStore.summarise_units(equation.lhs)
                    rhs_units = QuantityStore.summarise_units(equation.rhs)
                    assert QuantityStore.is_equal(lhs_units, rhs_units)

    def test_unit_extraction(self):
        eq = (5*units.mile/(2*units.hour + 10*units.minute))**(8*units.gram)
        assert QuantityStore.summarise_units(eq) == \
            (units.mile/(units.hour + units.minute))**units.gram

        millivolts = units.Quantity('millivolts', units.voltage, units.milli * units.volts, 'mV')
        x, y = sympy.symbols('x y')
        eq = (millivolts / units.millisecond)*sympy.Derivative(x, y)
        assert QuantityStore.summarise_units(eq) == (millivolts / units.milliseconds)

    def test_print(self, model):
        from os import environ
        if "CMLM_TEST_PRINT" in environ:
            # show equations
            for name, component in model.components.items():
                print(name)
                for equation in component.equations:
                    print('\t', equation)
