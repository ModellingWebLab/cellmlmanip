import os

import pytest
import sympy

from cellmlmanip import load_model, parser


class TestParser(object):
    @pytest.fixture(scope="class")
    def parser_instance(self):
        """Parses example CellML and returns model"""
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "test_simple_odes.cellml"
        )
        p = parser.Parser(example_cellml)
        return p

    @pytest.fixture(scope="class")
    def model(self, parser_instance):
        model = parser_instance.parse()
        return model

    def test_component_count(self, parser_instance, model):
        # grep -c '<component ' test_simple_odes.cellml
        assert len(parser_instance.components) == 21

    def test_group_relationships(self, parser_instance, model):
        assert parser_instance.components['circle_parent'].parent is None

        assert 'circle_x' in parser_instance.components['circle_parent'].encapsulated
        assert 'circle_y' in parser_instance.components['circle_parent'].encapsulated

        assert 'circle_x_source' in parser_instance.components['circle_x'].encapsulated
        assert 'circle_x_source' in parser_instance.components['circle_x_sibling'].siblings
        assert 'circle_x_sibling' in parser_instance.components['circle_x_source'].siblings
        assert 'circle_x' == parser_instance.components['circle_x_sibling'].parent

        assert 'circle_y_implementation' not in \
               parser_instance.components['circle_parent'].encapsulated

        assert 'circle_parent' == parser_instance.components['circle_x'].parent
        assert 'circle_parent' == parser_instance.components['circle_y'].parent

    def test_equations_count(self, model):
        # Include 19 equations from model + 2 equations added by parser for unit conversion
        assert len(model.equations) == 21  # NOTE: determined by eye!

    def test_variable_find(self, model):
        match = model.find_variable({'cmeta_id': 'time'})
        assert len(match) == 1
        assert match[0].cmeta_id == 'time' and match[0].name == 'environment$time'

        match = model.find_variable({'cmeta_id': 'sv12'})
        assert len(match) == 1 and \
            match[0].dummy.name == 'single_ode_rhs_const_var$sv1'

    def test_rdf(self, model):
        assert len(model.rdf) == 17

    def test_connections_loaded(self, model):
        # these are all unconnected variables
        unconnected = ['environment$time',
                       'single_independent_ode$sv1',
                       'single_ode_rhs_const_var$sv1',
                       'single_ode_rhs_const_var$a',
                       'single_ode_rhs_computed_var$sv1',
                       'single_ode_rhs_computed_var$a',
                       'derived_from_state_var$dbl_sv1',
                       'deriv_on_rhs$sv1_rate',
                       'circle_x_source$x',
                       'circle_x_sibling$x2',
                       'circle_y_implementation$y',
                       'circle_y_implementation$rhs',
                       'circle_sibling$local_complex_maths',
                       'time_units_conversion1$sv1',
                       'deriv_on_rhs1a$sv1_rate',
                       'time_units_conversion2$sv1',
                       'deriv_on_rhs2a$sv1_rate',
                       'state_units_conversion1$sv1',
                       'deriv_on_rhs1b$sv1_rate',
                       'state_units_conversion2$sv1',
                       'deriv_on_rhs2b$sv1_rate']
        for name in unconnected:
            variable = model.get_dummy_data(name)
            assert variable.dummy == variable.assigned_to

    def test_connections(self, model):
        # Check environment component's time variable has propagated
        environment__time = model.get_dummy_data('environment$time').dummy

        # We're checking sympy.Dummy objects (same name != same hash)
        assert isinstance(environment__time, sympy.Dummy)
        assert environment__time != sympy.Dummy(environment__time.name)

        state_units_conversion2__time = \
            model.get_dummy_data('state_units_conversion2$time').assigned_to
        assert environment__time == state_units_conversion2__time

        # propagated environment time to inside nested component circle_y
        circle_y__time = model.get_dummy_data('circle_y$time').assigned_to
        assert environment__time == circle_y__time

        # we have a new equation that links together times in different units
        time_units_conversion2__time = \
            model.get_dummy_data('time_units_conversion2$time').assigned_to
        equation = sympy.Eq(time_units_conversion2__time, environment__time)
        assert equation in model.equations

    def test_add_units_to_equations(self, model):
        # Eq(Derivative(_single_independent_ode$sv1, _environment$time), _1.00000000000000)
        # mV/millisecond == mV_per_ms
        test_equation = model.equations[0]
        lhs_units = model.units.summarise_units(test_equation.lhs)
        rhs_units = model.units.summarise_units(test_equation.rhs)
        assert model.units.is_unit_equal(rhs_units, lhs_units)

        # TODO: We should find two equations with different lhs/rhs units
        # 1. time_units_conversion1
        #    Eq(_time_units_conversion1$time, _environment$time) second millisecond
        # 2. time_units_conversion2
        #    Eq(_time_units_conversion2$time, _environment$time) microsecond millisecond
        # Make test to check two are unequal, fix them, then check equal

        # Try fixing all units on the RHS so that they match the LHS
        invalid_rhs_lhs_count = 0
        for index, equation in enumerate(model.equations):
            lhs_units = model.units.summarise_units(equation.lhs)
            rhs_units = model.units.summarise_units(equation.rhs)
            if not model.units.is_unit_equal(lhs_units, rhs_units):
                invalid_rhs_lhs_count += 1
                new_rhs = model.units.convert_to(1*rhs_units, lhs_units)
                # Create a new equality with the converted RHS and replace original
                conversion_factor = sympy.Dummy(str(new_rhs.magnitude))
                model.add_number(conversion_factor,
                                 {'cellml:units': str(lhs_units/rhs_units),
                                  'sympy.Number': sympy.Float(new_rhs.magnitude)}
                                 )
                equation = sympy.Eq(equation.lhs, equation.rhs * conversion_factor)
                # Replace the current equation with the same equation multiplied by factor
                model.equations[index] = equation
                lhs_units = model.units.summarise_units(equation.lhs)
                rhs_units = model.units.summarise_units(equation.rhs)
                # TODO: how to test this?
                assert model.units.is_unit_equal(lhs_units, rhs_units)
        assert invalid_rhs_lhs_count == 2

    @pytest.mark.skipif('CMLM_TEST_PRINT' not in os.environ, reason="print eq on demand")
    def test_print_eq(self, model):
        print()
        from cellmlmanip.units import ExpressionWithUnitPrinter
        symbol_info = dict()
        for var in model.dummy_data.values():
            symbol_info[var.dummy] = {'units': var.units}
            if var.number is not None:
                symbol_info[var.dummy]['number'] = var.number
        printer = ExpressionWithUnitPrinter(symbol_info=symbol_info)

        # show metadata for dummy instances in equations
        for index, key in enumerate(model.dummy_data.keys()):
            print('%3d. %s' % (index, str(model.dummy_data[key])))

        # show equations
        for index, equation in enumerate(model.equations):
            print('%3d. Eq(%s, %s)' % (index + 1,
                                       printer.doprint(equation.lhs),
                                       printer.doprint(equation.rhs)))
            lhs_units = model.units.summarise_units(equation.lhs)
            rhs_units = model.units.summarise_units(equation.rhs)
            print('     %s %s %s' %
                  (lhs_units,
                   '==' if model.units.is_unit_equal(rhs_units, lhs_units) else '!=',
                   rhs_units))

        dummy_instances, not_found = model.check_dummy_instances()
        assert len(not_found) == 0

    def test_connect_to_hidden_component(self):
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "err_connect_to_hidden_component.cellml"
        )
        p = parser.Parser(example_cellml)

        with pytest.raises(ValueError) as value_info:
            p.parse()

        assert 'Cannot determine the source & target' in str(value_info.value)

    def test_bad_connection_units(self):
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "err_bad_connection_units.cellml"
        )
        p = parser.Parser(example_cellml)
        model = p.parse()

        # then check the lhs/rhs units
        with pytest.raises(AssertionError) as assert_info:
            for e in model.equations:
                model.check_left_right_units_equal(e)

        match = ("Units volt (1.0, <Unit('kilogram * meter ** 2 / ampere / second ** 3')>) != "
                 "second (1.0, <Unit('second')>)")
        assert match in str(assert_info.value)

    def test_algebraic(self):
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "algebraic.cellml"
        )
        p = parser.Parser(example_cellml)
        model = p.parse()
        for e in model.equations:
            model.check_left_right_units_equal(e)

        ureg = model.units.ureg
        assert str(ureg.get_dimensionality(ureg.new_base)) == '[new_base]'
        assert ureg.get_base_units(ureg.new_base/ureg.second) == ureg.get_base_units(ureg.derived)

    def test_undefined_variable(self):
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "undefined_variable.cellml"
        )
        p = parser.Parser(example_cellml)
        with pytest.raises(AssertionError) as assert_info:
            p.parse()

        match = 'c$b not found in symbol dict'
        assert match in str(assert_info.value)

    def test_multiple_math_elements(self):
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "3.4.2.1.component_with_maths.cellml"
        )
        model = load_model(example_cellml)
        assert len(list(model.equations)) == 2
