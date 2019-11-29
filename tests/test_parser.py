import os

import pytest
import sympy

from cellmlmanip import load_model, parser

from . import shared


class TestParser(object):
    @pytest.fixture(scope="class")
    def parser_instance(self):
        """Returns an instance of the parser reading test_simple_odes.cellml"""
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "test_simple_odes.cellml"
        )
        p = parser.Parser(example_cellml)
        p.parse()
        return p

    def test_component_count(self, parser_instance):
        # grep -c '<component ' test_simple_odes.cellml
        assert len(parser_instance.components) == 21

    def test_group_relationships(self, parser_instance):
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

    def test_equations_count(self, simple_ode_model):
        # Include 19 equations from model + 2 equations added by parser for unit conversion
        assert len(simple_ode_model.equations) == 21  # NOTE: determined by eye!

    def test_rdf(self, simple_ode_model):
        assert len(simple_ode_model.rdf) == 21

    def test_connections_loaded(self, simple_ode_model):
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
            variable = simple_ode_model.get_symbol_by_name(name)
            assert variable == variable.assigned_to

    def test_connections(self, simple_ode_model):
        # Check environment component's time variable has propagated
        environment__time = simple_ode_model.get_symbol_by_name('environment$time')

        # We're checking sympy.Dummy objects (same name != same hash)
        assert isinstance(environment__time, sympy.Dummy)
        assert environment__time != sympy.Dummy(environment__time.name)

        state_units_conversion2__time = \
            simple_ode_model.get_symbol_by_name('state_units_conversion2$time').assigned_to
        assert environment__time == state_units_conversion2__time

        # propagated environment time to inside nested component circle_y
        circle_y__time = simple_ode_model.get_symbol_by_name('circle_y$time').assigned_to
        assert environment__time == circle_y__time

        # we have a new equation that links together times in different units
        time_units_conversion2__time = \
            simple_ode_model.get_symbol_by_name('time_units_conversion2$time').assigned_to
        # equation = sympy.Eq(time_units_conversion2__time, environment__time)
        equation = [e for e in simple_ode_model.equations if e.lhs == time_units_conversion2__time]
        assert len(equation) == 1
        equation = equation.pop()
        assert len(equation.rhs.find(environment__time)) == 1

    def test_add_units_to_equations(self, simple_ode_model):
        # Eq(Derivative(_single_independent_ode$sv1, _environment$time), _1.00000000000000)
        # mV/millisecond == mV_per_ms
        test_equation = simple_ode_model.equations[0]
        lhs_units = simple_ode_model.units.summarise_units(test_equation.lhs)
        rhs_units = simple_ode_model.units.summarise_units(test_equation.rhs)
        assert simple_ode_model.units.is_unit_equal(rhs_units, lhs_units)

        # two connected variables required conversion:
        # 1. time_units_conversion1
        #    Eq(_time_units_conversion1$time, _environment$time) second != millisecond
        # 2. time_units_conversion2
        #    Eq(_time_units_conversion2$time, _environment$time) microsecond != millisecond

        require_conversion = [simple_ode_model.get_symbol_by_name('time_units_conversion1$time').assigned_to,
                              simple_ode_model.get_symbol_by_name('time_units_conversion2$time').assigned_to]

        # find the equations that define these variables that require conversion
        invalid_rhs_lhs_count = 0
        for index, equation in enumerate(simple_ode_model.equations):
            if equation.lhs in require_conversion:
                # the lhs and rhs units should be equal
                invalid_rhs_lhs_count += 1
                lhs_units = simple_ode_model.units.summarise_units(equation.lhs)
                rhs_units = simple_ode_model.units.summarise_units(equation.rhs)
                assert simple_ode_model.units.is_unit_equal(lhs_units, rhs_units)
        assert invalid_rhs_lhs_count == 2

    @pytest.mark.skipif('CMLM_TEST_PRINT' not in os.environ, reason="print eq on demand")
    def test_print_eq(self, simple_ode_model):
        print()
        from cellmlmanip.units import ExpressionWithUnitPrinter
        printer = ExpressionWithUnitPrinter()

        # show equations
        for index, equation in enumerate(simple_ode_model.equations):
            print('%3d. Eq(%s, %s)' % (index + 1,
                                       printer.doprint(equation.lhs),
                                       printer.doprint(equation.rhs)))
            lhs_units = simple_ode_model.units.summarise_units(equation.lhs)
            rhs_units = simple_ode_model.units.summarise_units(equation.rhs)
            print('     %s %s %s' %
                  (lhs_units,
                   '==' if simple_ode_model.units.is_unit_equal(rhs_units, lhs_units) else '!=',
                   rhs_units))

        assert shared.check_cmeta_ids(simple_ode_model)
        assert shared.check_dummy_assignment(simple_ode_model)

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

        from pint.errors import DimensionalityError
        with pytest.raises(DimensionalityError) as dim_error:
            p.parse()

        match = ("Cannot convert from 'second' ([time]) to "
                 "'volt' ([length] ** 2 * [mass] / [current] / [time] ** 3)")

        assert match in str(dim_error.value)

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
        assert ureg.get_base_units(ureg.new_base / ureg.second) == ureg.get_base_units(ureg.derived)

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

    def test_bad_unit_definitions(self):
        bad_unit_models = ['5.4.2.2.unit_cycle_1.cellml',
                           '5.4.2.2.unit_cycle_3.cellml',
                           '5.4.2.2.unit_units_invalid.cellml']
        for bad_unit_model in bad_unit_models:
            example_cellml = os.path.join(
                os.path.dirname(__file__), "cellml_files", bad_unit_model
            )
            p = parser.Parser(example_cellml)
            with pytest.raises(ValueError) as value_error:
                p.parse()
            assert 'Cannot create units' in str(value_error.value)

    def test_component_units_unsupported(self):
        """ Component units are not supported. """

        path = os.path.join(os.path.dirname(__file__), 'cellml_files', '3.4.2.1.component_with_units.cellml')
        p = parser.Parser(path)
        with pytest.raises(ValueError, match='Defining units inside components is not supported'):
            p.parse()

    def test_reactions_unsupported(self):
        """ Reactions are not supported. """

        path = os.path.join(os.path.dirname(__file__), 'cellml_files', '3.4.2.1.component_with_reactions.cellml')
        p = parser.Parser(path)
        with pytest.raises(ValueError, match='Reactions are not supported'):
            p.parse()

