import os

import pytest
import sympy

from cellmlmanip import load_model, parser
from cellmlmanip.model import NumberWrapper


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
        print(len(model.equations_x))
        # Include 19 equations from model + 2 equations added by parser for unit conversion
        assert len(model.equations_x) == 21  # NOTE: determined by eye!

    def test_variable_find(self, model):
        match = model.find_variable_x({'cmeta_id': 'time'})
        assert len(match) == 1
        assert match[0].cmeta_id == 'time' and match[0].name == 'environment$time'

        match = model.find_variable_x({'cmeta_id': 'sv12'})
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
            variable = model.variables_x[name]
            assert variable.dummy == variable.assignment

    def test_connections(self, model):
        # Check environment component's time variable has propagated
        environment__time = model.variables_x['environment$time'].dummy

        # We're checking sympy.Dummy objects (same name != same hash)
        assert isinstance(environment__time, sympy.Dummy)
        assert environment__time != sympy.Dummy(environment__time.name)

        state_units_conversion2__time = \
            model.variables_x['state_units_conversion2$time'].assignment
        assert environment__time == state_units_conversion2__time

        # propagated environment time to inside nested component circle_y
        circle_y__time = model.variables_x['circle_y$time'].assignment
        assert environment__time == circle_y__time

        # we have a new equation that links together times in different units
        time_units_conversion2__time = \
            model.variables_x['time_units_conversion2$time'].assignment
        equation = sympy.Eq(time_units_conversion2__time, environment__time)
        assert equation in model.equations_x

    def test_add_units_to_equations(self, model):
        # Eq(Derivative(_single_independent_ode$sv1, _environment$time), _1.00000000000000)
        # mV/millisecond == mV_per_ms
        test_equation = model.equations_x[0]
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
        for index, equation in enumerate(model.equations_x):
            lhs_units = model.units.summarise_units(equation.lhs)
            rhs_units = model.units.summarise_units(equation.rhs)
            if not model.units.is_unit_equal(lhs_units, rhs_units):
                invalid_rhs_lhs_count += 1
                new_rhs = model.units.convert_to(1*rhs_units, lhs_units)
                # Create a new equality with the converted RHS and replace original
                new_dummy = sympy.Dummy(str(new_rhs.magnitude))
                model.numbers_x[new_dummy] = NumberWrapper(((1*lhs_units) / (1*rhs_units)).units,
                                                           sympy.Float(new_rhs.magnitude))
                equation = sympy.Eq(equation.lhs, equation.rhs * new_dummy)
                # Replace the current equation with the same equation multiplied by factor
                model.equations_x[index] = equation
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
        for var in model.variables_x.values():
            symbol_info[var.dummy] = {'units': var.units}
        for d, num in model.numbers_x.items():
            symbol_info[d] = {'units': num.units, 'number': num.number}
        printer = ExpressionWithUnitPrinter(symbol_info=symbol_info)
        # show equations
        for index, equation in enumerate(model.equations_x):
            print('%3d. Eq(%s, %s)' % (index + 1,
                                       printer.doprint(equation.lhs),
                                       printer.doprint(equation.rhs)))
            lhs_units = model.units.summarise_units(equation.lhs)
            rhs_units = model.units.summarise_units(equation.rhs)
            print('     %s %s %s' %
                  (lhs_units,
                   '==' if model.units.is_unit_equal(rhs_units, lhs_units) else '!=',
                   rhs_units))

    def test_connect_to_hidden_component(self):
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "err_connect_to_hidden_component.cellml"
        )
        p = parser.Parser(example_cellml)

        with pytest.raises(ValueError) as value_info:
            model = p.parse()
            model.make_connections()

        assert 'Cannot determine the source & target' in str(value_info.value)

    def test_bad_connection_units(self):
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "err_bad_connection_units.cellml"
        )
        p = parser.Parser(example_cellml)
        model = p.parse()

        # then check the lhs/rhs units
        with pytest.raises(AssertionError) as assert_info:
            for e in model.equations_x:
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
        for e in model.equations_x:
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
        assert len(list(model.equations_x)) == 2
