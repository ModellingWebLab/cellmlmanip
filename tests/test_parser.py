import logging
import os

import pytest
import sympy

from cellmlmanip import parser

from .shared import check_left_right_units_equal, load_model


logger = logging.getLogger(__name__)


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
        """ Tests correct number of components parsed."""
        # grep -c '<component ' test_simple_odes.cellml
        assert len(parser_instance.components) == 21

    def test_group_relationships(self, parser_instance):
        """ Tests the correct relationships are parsed."""
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
        """Tests correct number of equations read and then created where necessary by parser. """
        # Include 19 equations from model + 2 equations added by parser for unit conversion + 1 constant
        assert len(simple_ode_model.equations) == 22  # NOTE: determined by eye!

    def test_rdf(self, simple_ode_model):
        """ Tests correct number of rdf statements parsed. """
        assert len(simple_ode_model.rdf) == 21

    def test_connections_loaded(self, simple_ode_model):
        """ Tests the variables that are assigned to themselves.
        Note these variables do not occur on the lhs of any equation."""
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
            variable = simple_ode_model.get_variable_by_name(name)
            assert variable == variable.assigned_to

    def test_connections(self, simple_ode_model):
        """ Tests the variables that are assigned to other variable.
        Note these variables appear on lhs of an equation.
        """
        # Check environment component's time variable has propagated
        environment__time = simple_ode_model.get_variable_by_name('environment$time')

        # We're checking sympy.Dummy objects (same name != same hash)
        assert isinstance(environment__time, sympy.Dummy)
        assert environment__time != sympy.Dummy(environment__time.name)

        state_units_conversion2__time = \
            simple_ode_model.get_variable_by_name('state_units_conversion2$time').assigned_to
        assert environment__time == state_units_conversion2__time

        # propagated environment time to inside nested component circle_y
        circle_y__time = simple_ode_model.get_variable_by_name('circle_y$time').assigned_to
        assert environment__time == circle_y__time

        # we have a new equation that links together times in different units
        time_units_conversion2__time = \
            simple_ode_model.get_variable_by_name('time_units_conversion2$time').assigned_to
        # equation = sympy.Eq(time_units_conversion2__time, environment__time)
        equation = [e for e in simple_ode_model.equations if e.lhs == time_units_conversion2__time]
        assert len(equation) == 1
        equation = equation.pop()
        assert len(equation.rhs.find(environment__time)) == 1

    def test_add_units_to_equations(self, simple_ode_model):
        """ Tests that parser reads all equations, adds and works out units and then
        adds equations for necessary conversions"""
        # Eq(Derivative(_single_independent_ode$sv1, _environment$time), _1.00000000000000)
        # mV/millisecond == mV_per_ms
        test_equation = simple_ode_model.equations[0]
        lhs_units = simple_ode_model.units.evaluate_units(test_equation.lhs)
        rhs_units = simple_ode_model.units.evaluate_units(test_equation.rhs)
        assert simple_ode_model.units.is_equivalent(rhs_units, lhs_units)

        # two connected variables required conversion:
        # 1. time_units_conversion1
        #    Eq(_time_units_conversion1$time, _environment$time) second != millisecond
        # 2. time_units_conversion2
        #    Eq(_time_units_conversion2$time, _environment$time) microsecond != millisecond

        require_conversion = [simple_ode_model.get_variable_by_name('time_units_conversion1$time').assigned_to,
                              simple_ode_model.get_variable_by_name('time_units_conversion2$time').assigned_to]

        # find the equations that define these variables that require conversion
        invalid_rhs_lhs_count = 0
        for index, equation in enumerate(simple_ode_model.equations):
            if equation.lhs in require_conversion:
                # the lhs and rhs units should be equal
                invalid_rhs_lhs_count += 1
                lhs_units = simple_ode_model.units.evaluate_units(equation.lhs)
                rhs_units = simple_ode_model.units.evaluate_units(equation.rhs)
                assert simple_ode_model.units.is_equivalent(lhs_units, rhs_units)
        assert invalid_rhs_lhs_count == 2

    def test_check_rdf_identities(self, simple_ode_model):
        """Checks that every variable in a model with an ``rdf_identity`` is a source variable."""
        for variable in simple_ode_model.graph:
            if variable.is_Derivative:
                variable = variable.free_symbols.pop()
            if variable.rdf_identity is not None:
                assert variable == variable.assigned_to, ('%s has cmeta id but is assigned to %s',
                                                          variable.dummy, variable.assigned_to)

    def check_dummy_assignment(self, simple_ode_model):
        """Every variable in the model should be assigned to itself or another variable.
        The source variable must be assigned to itself."""
        for variable in simple_ode_model.graph:
            if variable.is_Derivative:
                variable = variable.free_symbols.pop()
            # either the variable is assigned to itself
            if variable == variable.assigned_to:
                continue

            # or the variable is assigned to a source variable
            source = variable.assigned_to

            # the source dummy must be assigned to itself
            assert source.assigned_to == source, ('%s is assigned to %s, which is assigned to %s',
                                                  variable, source, source.assigned_to)

    def test_connect_to_hidden_component(self):
        """ Tests parser throws exception when cannot connect components. """
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "err_connect_to_hidden_component.cellml"
        )
        p = parser.Parser(example_cellml)

        with pytest.raises(ValueError) as value_info:
            p.parse()

        assert 'Cannot determine the source & target' in str(value_info.value)

    def test_bad_connection_units(self):
        """ Tests parser throws exception when a units of connected variables
         are incompatible and a conversion equation cannot be created. """
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

    def test_new_base_units(self):
        """Tests unit checking in a model that defines a new base unit."""
        cellml = os.path.join(os.path.dirname(__file__), "cellml_files", "algebraic.cellml")
        p = parser.Parser(cellml)
        model = p.parse()
        for e in model.equations:
            check_left_right_units_equal(model.units, e)

    def test_algebraic_model(self):
        """Tests parsing of a model with no ODEs works."""
        model = load_model('algebraic')

        # Getting the free variable raises an exception
        with pytest.raises(ValueError):
            model.get_free_variable()

    def test_units_with_multiplier(self):
        """Tests parsing a unit with a multiplier."""
        cellml = os.path.join(os.path.dirname(__file__), 'cellml_files', 'imperial_units.cellml')
        p = parser.Parser(cellml)
        m = p.parse()
        i = m.units.get_unit('inch')
        assert m.units.format(i, True) == '0.0254 meter'

    def test_undefined_variable(self):
        """ Tests parser exception for undefined variable. """
        example_cellml = os.path.join(
            os.path.dirname(__file__), "cellml_files", "undefined_variable.cellml"
        )
        p = parser.Parser(example_cellml)
        with pytest.raises(AssertionError) as assert_info:
            p.parse()

        match = 'c$b not found in symbol dict'
        assert match in str(assert_info.value)

    def test_multiple_math_elements(self):
        """ Tests parser correctly handles a component with more than one math element. """
        model = load_model('3.4.2.1.component_with_maths.cellml')
        assert len(list(model.equations)) == 2

    def test_bad_unit_definitions(self):
        """ Tests parser exception for invalid unit declarations in a model. """
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
        """ Tests parser exception: Component units are not supported. """

        path = os.path.join(os.path.dirname(__file__), 'cellml_files', '3.4.2.1.component_with_units.cellml')
        p = parser.Parser(path)
        with pytest.raises(ValueError, match='Defining units inside components is not supported'):
            p.parse()

    def test_reactions_unsupported(self):
        """ Tests parser exception: Reactions are not supported. """

        path = os.path.join(os.path.dirname(__file__), 'cellml_files', '3.4.2.1.component_with_reactions.cellml')
        p = parser.Parser(path)
        with pytest.raises(ValueError, match='Reactions are not supported'):
            p.parse()

    def test_transform_constants(self):
        """ Tests Parser.transform_constants(). """

        # Parse model, which should convert an initial value to an equation
        model = load_model('initial_value_constant.cellml')
        v = model.get_variable_by_name('A$a')

        # Initial value should be None, if transform_constants worked
        assert v.initial_value is None

        # And its RHS value should be 1
        assert model.get_value(v) == 1

        # And the equations should have matching units
        for eq in model.equations:
            check_left_right_units_equal(model.units, eq)

    def test_validation(self):
        """ Tests that RNG validation can pick stuff up. """

        path = os.path.join(os.path.dirname(__file__), 'cellml_files', 'err_extra_content.cellml')
        p = parser.Parser(path)
        with pytest.raises(ValueError, match='Element model has extra content'):
            p.parse()

    def test_bad_encapsulation(self):
        with pytest.raises(ValueError, match="Encapsulated component circle_x already added!"):
            load_model('test_simple_odes_bad_encapsulation')

    def test_wrong_parent_component(self):
        with pytest.raises(ValueError, match=("Parent of component dummy2 already dummy. "
                                              "Cannot set environment: multiple parents not allowed!")):
            load_model('test_wrong_parent_component')

    def test_wrong_encapsulation(self):
        with pytest.raises(ValueError, match="Encapsulated component circle_x_sibling already added!"):
            load_model('test_wrong_encapsulation')

    def test_wrong_relationship_tags(self):
        with pytest.raises(ValueError, match="Expecting exactly 1 relationship_ref tag per group, got 2!"):
            load_model('test_wrong_relationship_tags')

    def test_wrong_connections(self):
        with pytest.raises(ValueError, match=(r'Target already assigned to environment\$time '
                                              r'before assignment to environment\$time')):
            load_model('test_wrong_connections')

    def test_piecewise_booleans_error(self):
        # Testing a model with a boolean in expression being passed to Piecewise
        # see https://github.com/sympy/sympy/issues/24086
        # and https://github.com/ModellingWebLab/cellmlmanip/issues/350
        model = load_model('parsing_err_bool_in_cond.cellml')
        assert sorted(map(str, model.variables())) == ['A$iffalse', 'A$iftrue', 'A$x']

        zero, one, Ax, iffalse, iftrue = sympy.symbols('_0, _1, _A$x, _A$iffalse, _A$iftrue')
        eqs1 = [sympy.Eq(iffalse, sympy.Piecewise((one, sympy.Eq(False, sympy.Eq(zero, Ax))), (zero, True))),
                sympy.Eq(iftrue, sympy.Piecewise((one, sympy.Eq(True, sympy.Eq(zero, Ax))), (zero, True))),
                sympy.Eq(Ax, zero)]
        eqs2 = [sympy.Eq(iffalse, sympy.Piecewise((one, sympy.Eq(False, sympy.Eq(Ax, zero))), (zero, True))),
                sympy.Eq(iftrue, sympy.Piecewise((one, sympy.Eq(True, sympy.Eq(Ax, zero))), (zero, True))),
                sympy.Eq(Ax, zero)]
        model_eqs = sorted(map(str, model.equations))
        assert model_eqs == sorted(map(str, set(eqs1))) or model_eqs == sorted(map(str, set(eqs2)))
