import logging
import os

import pytest
import sympy

from cellmlmanip import parser

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

    def test_group_relationships(self, model):
        assert model.components['circle_parent'].parent is None

        assert 'circle_x' in model.components['circle_parent'].encapsulated
        assert 'circle_y' in model.components['circle_parent'].encapsulated

        assert 'circle_parent' == model.components['circle_x'].parent
        assert 'circle_parent' == model.components['circle_y'].parent

        assert 'circle_x_source' in model.components['circle_x'].encapsulated
        assert 'circle_x_source' in model.components['circle_x_sibling'].siblings
        assert 'circle_x_sibling' in model.components['circle_x_source'].siblings
        assert 'circle_x' == model.components['circle_x_sibling'].parent

        assert 'circle_y_implementation' not in model.components['circle_parent'].encapsulated

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
            matched[0]['sympy.Dummy'].name == 'single_ode_rhs_const_var$sv1'

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

    def test_add_units_to_equations(self, model):
        # This is an irreversible operation # TODO: don't mutate?
        model.add_units_to_equations()

        # mV/millisecond == mV_per_ms
        test_equation = model.components['single_independent_ode'].equations[0]
        lhs_units = model.units.summarise_units(test_equation.lhs)
        rhs_units = model.units.summarise_units(test_equation.rhs)
        assert model.units.is_unit_equal(rhs_units, lhs_units)

        # mV_per_usec != uV/millisecond
        test_equation = model.components['deriv_on_rhs2b'].equations[0]
        lhs_units = model.units.summarise_units(test_equation.lhs)
        rhs_units = model.units.summarise_units(test_equation.rhs)
        assert not model.units.is_unit_equal(rhs_units, lhs_units)

        # Convert RHS into required units. This returns a Quantity object
        new_rhs = model.units.convert_to(rhs_units, lhs_units)
        assert new_rhs.magnitude == pytest.approx(1.0e-6)
        assert model.units.is_unit_equal(new_rhs.units, lhs_units)

        # TODO: try conversion of units of RHS by LHS.units / RHS.unit == x;
        # i.e. if x == 1, good, else RHS = RHS * x

        # Try fixing all units on the RHS so that they match the LHS
        for component in model.components.values():
            for index, equation in enumerate(component.equations):
                lhs_units = model.units.summarise_units(equation.lhs)
                rhs_units = model.units.summarise_units(equation.rhs)
                if not model.units.is_unit_equal(lhs_units, rhs_units):
                    new_rhs = model.units.convert_to(rhs_units, lhs_units)
                    # Create a new equality with the converted RHS and replace original
                    from cellmlmanip.units import UnitDummy
                    new_dummy = UnitDummy(str(new_rhs.magnitude))
                    new_dummy.number = sympy.Float(new_rhs.magnitude)
                    new_dummy.unit = ((1*lhs_units) / (1*rhs_units)).units
                    equation = sympy.Eq(equation.lhs, equation.rhs * new_dummy)
                    component.equations[index] = equation
                    print('new eq:', equation)
                    lhs_units = model.units.summarise_units(equation.lhs)
                    rhs_units = model.units.summarise_units(equation.rhs)
                    # TODO: how to test this?
                    assert model.units.is_unit_equal(lhs_units, rhs_units)

    @pytest.mark.skipif('CMLM_TEST_PRINT' not in os.environ, reason="print eq on demand")
    def test_print_eq(self, model):
        # show equations
        for name, component in model.components.items():
            print(name)
            for equation in component.equations:
                print('\t', equation)

    def test_connect_to_hidden_component(self):
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "err_connect_to_hidden_component.cellml"
        )
        p = parser.Parser(example_cellml)
        model = p.parse()
        with pytest.raises(ValueError, match=r'^Cannot determine the source & target.*'):
            model.make_connections()

    def test_bad_connection_units(self):
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "err_bad_connection_units.cellml"
        )
        p = parser.Parser(example_cellml)
        model = p.parse()

        # first we make the connections
        model.make_connections()

        # then add the units to the equations
        model.add_units_to_equations()

        # then check the lhs/rhs units
        with pytest.raises(AssertionError, match='Units second != volt'):
            for e in model.equations:
                model.check_left_right_units_equal(e)

    # def test_algebraic(self):
    #     example_cellml = os.path.join(
    #         os.path.dirname(__file__), "cellml_files", "algebraic.cellml"
    #     )
    #     p = parser.Parser(example_cellml)
    #     model = p.parse()
    #     model.make_connections()
    #     model.add_units_to_equations()
    #     for e in model.equations:
    #         model.check_left_right_units_equal(e)

    def test_undefined_variable(self):
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "undefined_variable.cellml"
        )
        p = parser.Parser(example_cellml)
        with pytest.raises(KeyError, match=r'Variable "b" in component "c" could not be found\.'):
            p.parse()
