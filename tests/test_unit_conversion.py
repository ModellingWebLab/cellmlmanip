import pytest
from pint import DimensionalityError

from cellmlmanip import units
from cellmlmanip.model import DataDirectionFlow
from . import shared


class TestUnitConversion:
    ###############################################################
    # fixtures

    @pytest.fixture
    def local_model(scope='function'):
        """ Fixture to load a local copy of  the basic_ode model that may get modified. """
        return shared.load_model('basic_ode')

    @pytest.fixture
    def br_model(scope='function'):
        """ Fixture to load a local copy of  the basic_ode model that may get modified. """
        return shared.load_model('beeler_reuter_model_1977')

    @pytest.fixture
    def literals_model(scope='function'):
        """ Fixture to load a local copy of  the basic_ode model that may get modified. """
        return shared.load_model('literals_for_conversion_tests')

    @pytest.fixture
    def multiode_model(scope='function'):
        """ Fixture to load a local copy of  the basic_ode model that may get modified. """
        return shared.load_model('repeated_ode_for_conversion_tests')

    @pytest.fixture
    def multiode_freevar_model(scope='function'):
        """ Fixture to load a local copy of  the basic_ode model that may get modified. """
        return shared.load_model('repeated_ode_freevar_for_conversion_tests')

    @pytest.fixture
    def silly_names(scope='function'):
        """ Fixture to load a local copy of  the basic_ode model that may get modified. """
        return shared.load_model('silly_names')

    ###############################################################
    # Basic tests

    def test_bad_units(self, bad_units_model):
        """ Tests units read and calculated from an inconsistent model. """
        symbol_a = bad_units_model.get_variable_by_cmeta_id("a")
        symbol_b = bad_units_model.get_variable_by_cmeta_id("b")
        equation = bad_units_model.get_equations_for([symbol_b], strip_units=False)
        assert len(equation) == 2
        assert equation[0].lhs == symbol_a
        assert bad_units_model.units.evaluate_units(equation[0].lhs) == bad_units_model.units.get_unit('ms')
        with pytest.raises(units.UnitError):
            # cellml file states a (ms) = 1 (ms) + 1 (second)
            bad_units_model.units.evaluate_units(equation[0].rhs)

        assert equation[1].lhs == symbol_b
        with pytest.raises(units.UnitError):
            # cellml file states b (per_ms) = power(a (ms), 1 (second))
            bad_units_model.units.evaluate_units(equation[1].rhs)

    def test_providing_unit_store(self):
        """Tests that we can provide an existing UnitStore to a model."""
        store = units.UnitStore()
        model = shared.load_model('basic_ode', unit_store=store)
        original_var = model.get_variable_by_name('env_ode$sv1')

        # test no change in units
        store.add_unit('my_mV', '1e-3 * volt')
        my_mV = store.get_unit('my_mV')
        new_var = model.convert_variable(original_var, my_mV, DataDirectionFlow.INPUT)
        assert new_var == original_var

        # test conversion
        store.add_unit('my_nV', '1e-9 * volt')
        my_nV = store.get_unit('my_nV')
        new_var = model.convert_variable(original_var, my_nV, DataDirectionFlow.OUTPUT)
        assert new_var != original_var
        assert str(model.get_definition(new_var)) == 'Eq(_env_ode$sv1_converted, 1000000.0*_env_ode$sv1)'

    ###############################################################
    # Helper functions for later tests

    # original state for local_model
    def _original_state_local_model(self, local_model):
        assert len(local_model.variables()) == 3
        symbol_a = local_model.get_variable_by_cmeta_id('sv11')
        symbol_t = local_model.get_variable_by_cmeta_id('time')
        assert symbol_a.initial_value == 2.0
        assert symbol_a.units == local_model.units.get_unit('mV')
        assert symbol_t.units == local_model.units.get_unit('ms')
        assert len(local_model.equations) == 1
        assert str(local_model.equations[0]) == 'Eq(Derivative(_env_ode$sv1, _environment$time), _1.0)'
        states = local_model.get_state_variables()
        assert len(states) == 1
        assert symbol_a in states
        return True

    # original state for literals_model
    def _original_state_literals_model(self, literals_model):
        assert len(literals_model.variables()) == 5
        symbol_a = literals_model.get_variable_by_cmeta_id('sv11')
        symbol_t = literals_model.get_variable_by_cmeta_id('time')
        symbol_x = literals_model.get_variable_by_cmeta_id('current')
        symbol_y = literals_model.get_variable_by_name('env_ode$y')
        assert symbol_x.name == 'env_ode$x'
        assert symbol_a.initial_value == 2.0
        assert symbol_t.initial_value is None
        assert symbol_x.initial_value is None
        assert symbol_y.initial_value is None
        assert symbol_a.units == literals_model.units.get_unit('mV')
        assert symbol_t.units == literals_model.units.get_unit('ms')
        assert symbol_x.units == literals_model.units.get_unit('pA')
        assert symbol_y.units == literals_model.units.get_unit('per_pA')
        assert len(literals_model.equations) == 3
        assert str(literals_model.equations[0]) == 'Eq(Derivative(_env_ode$sv1, _environment$time), _1.0)'
        assert str(literals_model.equations[1]) == 'Eq(_env_ode$x, _1.0)'
        assert str(literals_model.equations[2]) == 'Eq(_env_ode$y, _1.0/_env_ode$x)'
        return True

    # original state for multiode_model
    def _original_state_multiode_model(self, multiode_model):
        assert len(multiode_model.variables()) == 5
        symbol_a = multiode_model.get_variable_by_cmeta_id('sv11')
        symbol_t = multiode_model.get_variable_by_cmeta_id('time')
        assert symbol_a.initial_value == 2.0
        assert symbol_a.units == multiode_model.units.get_unit('mV')
        assert symbol_t.units == multiode_model.units.get_unit('ms')
        assert len(multiode_model.equations) == 3
        assert str(multiode_model.equations[0]) == 'Eq(Derivative(_env_ode$sv1, _environment$time), _1.0)'
        assert str(multiode_model.equations[1]) == \
            'Eq(_env_ode$x, _1.0 + _3.0*Derivative(_env_ode$sv1, _environment$time))'
        assert str(multiode_model.equations[2]) == 'Eq(Derivative(_env_ode$y, _environment$time), _2.0)'
        return True

    # original state
    def _original_state_multiode_freevar(self, multiode_freevar_model):
        assert len(multiode_freevar_model.variables()) == 4
        symbol_a = multiode_freevar_model.get_variable_by_cmeta_id('sv11')
        symbol_t = multiode_freevar_model.get_variable_by_cmeta_id('time')
        symbol_y = multiode_freevar_model.get_variable_by_name('env_ode$y')
        assert symbol_a.initial_value == 2.0
        assert symbol_y.initial_value == 3.0
        assert symbol_a.units == multiode_freevar_model.units.get_unit('mV')
        assert symbol_t.units == multiode_freevar_model.units.get_unit('ms')
        assert symbol_y.units == multiode_freevar_model.units.get_unit('mV')
        assert len(multiode_freevar_model.equations) == 2
        assert str(multiode_freevar_model.equations[0]) == 'Eq(Derivative(_env_ode$sv1, _environment$time), _1.0)'
        assert str(multiode_freevar_model.equations[1]) == 'Eq(Derivative(_env_ode$y, _environment$time), _2.0)'
        return True

    ###############################################################
    # Main conversion tests

    def test_add_input_state_variable(self, local_model):
        """ Tests the Model.convert_variable function that changes units.
        This particular test is when a state variable is being converted

        For example::

            Original model
                var{time} time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};

                ode(sv1, time) = 1{mV_per_ms};

        convert_variable(sv1, volt, DataDirectionFlow.INPUT)

            creates model
                var{time} time: ms {pub: in};
                var{sv11} sv1_converted: V {init: 0.002};
                var sv1 mV
                var sv1_orig_deriv mV_per_ms

                sv1 = 1000 * sv1_converted
                sv1_orig_deriv = 1{mV_per_ms}
                ode(sv1_converted, time) = 0.001 * sv1_orig_deriv
        """
        mV_unit = local_model.units.get_unit('mV')
        volt_unit = local_model.units.get_unit('volt')
        original_var = local_model.get_variable_by_name('env_ode$sv1')

        assert self._original_state_local_model(local_model)
        # test no change in units
        newvar = local_model.convert_variable(original_var, mV_unit, DataDirectionFlow.INPUT)
        assert newvar == original_var
        assert self._original_state_local_model(local_model)

        # change mV to V
        newvar = local_model.convert_variable(original_var, volt_unit, DataDirectionFlow.INPUT)
        assert newvar != original_var
        assert len(local_model.variables()) == 5
        symbol_a = local_model.get_variable_by_cmeta_id('sv11')
        assert symbol_a.initial_value == 0.002
        assert symbol_a.units == local_model.units.get_unit('volt')
        assert symbol_a.name == 'env_ode$sv1_converted'
        symbol_t = local_model.get_variable_by_cmeta_id('time')
        assert symbol_t.units == local_model.units.get_unit('ms')
        symbol_orig = local_model.get_variable_by_name('env_ode$sv1')
        assert symbol_orig.units == local_model.units.get_unit('mV')
        assert symbol_orig.initial_value is None
        symbol_derv = local_model.get_variable_by_name('env_ode$sv1_orig_deriv')
        assert symbol_derv.units == local_model.units.get_unit('mV') / local_model.units.get_unit('ms')
        assert symbol_derv.initial_value is None
        assert len(local_model.equations) == 3
        assert str(local_model.equations[0]) == 'Eq(_env_ode$sv1, 1000.0*_env_ode$sv1_converted)'
        assert str(local_model.equations[1]) == 'Eq(_env_ode$sv1_orig_deriv, _1.0)'
        assert str(local_model.equations[2]) == 'Eq(Derivative(_env_ode$sv1_converted, _environment$time), ' \
                                                '0.001*_env_ode$sv1_orig_deriv)'

        states = local_model.get_state_variables()
        assert len(states) == 1
        assert symbol_a in states

    def test_add_input_free_variable(self, local_model):
        """ Tests the Model.convert_variable function that changes units of given variable.
        This particular case tests changing a free variable
        For example::

            Original model
                var{time} time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};

                ode(sv1, time) = 1{mV_per_ms};

        convert_variable(time, second, DataDirectionFlow.INPUT)

        becomes
                var time: ms;
                var{time} time_converted: s;
                var{sv11} sv1: mV {init: 2};
                var sv1_orig_deriv mV_per_ms

                time = 1000 * time_converted;
                sv1_orig_deriv = 1{mV_per_ms}
                ode(sv1, time_converted) = 1000 * sv1_orig_deriv
        """
        ms_unit = local_model.units.get_unit('ms')
        second_unit = local_model.units.get_unit('second')
        original_var = local_model.get_variable_by_name('environment$time')

        assert self._original_state_local_model(local_model)
        # test no change in units
        local_model.convert_variable(original_var, ms_unit, DataDirectionFlow.INPUT)
        assert self._original_state_local_model(local_model)

        # change ms to s
        local_model.convert_variable(original_var, second_unit, DataDirectionFlow.INPUT)
        assert len(local_model.variables()) == 5
        symbol_a = local_model.get_variable_by_cmeta_id('sv11')
        assert symbol_a.initial_value == 2.0
        assert symbol_a.units == local_model.units.get_unit('mV')
        assert symbol_a.name == 'env_ode$sv1'
        symbol_t = local_model.get_variable_by_cmeta_id('time')
        assert symbol_t.units == local_model.units.get_unit('second')
        assert symbol_t.name == 'environment$time_converted'
        symbol_orig = local_model.get_variable_by_name('env_ode$sv1')
        assert symbol_orig.units == local_model.units.get_unit('mV')
        symbol_derv = local_model.get_variable_by_name('env_ode$sv1_orig_deriv')
        assert symbol_derv.units == local_model.units.get_unit('mV') / local_model.units.get_unit('ms')
        assert len(local_model.equations) == 3
        assert str(local_model.equations[0]) == 'Eq(_environment$time, 1000.0*_environment$time_converted)'
        assert str(local_model.equations[1]) == 'Eq(_env_ode$sv1_orig_deriv, _1.0)'
        assert str(local_model.equations[2]) == 'Eq(Derivative(_env_ode$sv1, _environment$time_converted), ' \
                                                '1000.0*_env_ode$sv1_orig_deriv)'

    def test_add_input_literal_variable(self, literals_model):
        """ Tests the Model.convert_variable function that changes units of given variable.
        This particular case tests changing a literal variable/constant
        For example::

            Original model
                var{time} time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};
                var{current} x: pA;
                var y: per_pA;

                ode(sv1, time) = 1{mV_per_ms};
                x = 1{pA};
                y = 1{dimensionless}/x;

        convert_variable(x, nA, DataDirectionFlow.INPUT)

        becomes
                var time: ms;
                var{sv11} sv1: mV {init: 2};
                var{current} x_converted: nA;
                var x : pA
                var y: per_pA;

                ode(sv1, time) = 1{mV_per_ms};
                x_converted = 0.001 * 1{pA}
                x = 1000* x_converted;
                y = 1{dimensionless}/x;
        """
        pA_unit = literals_model.units.get_unit('pA')
        nA_unit = literals_model.units.get_unit('nA')
        original_var = literals_model.get_variable_by_name('env_ode$x')

        assert self._original_state_literals_model(literals_model)
        # test no change in units
        literals_model.convert_variable(original_var, pA_unit, DataDirectionFlow.INPUT)
        assert self._original_state_literals_model(literals_model)

        # change pA to nA
        literals_model.convert_variable(original_var, nA_unit, DataDirectionFlow.INPUT)
        assert len(literals_model.variables()) == 6
        symbol_a = literals_model.get_variable_by_cmeta_id('sv11')
        symbol_t = literals_model.get_variable_by_cmeta_id('time')
        symbol_x = literals_model.get_variable_by_cmeta_id('current')
        assert symbol_x.name == 'env_ode$x_converted'
        symbol_y = literals_model.get_variable_by_name('env_ode$y')
        symbol_x_orig = literals_model.get_variable_by_name('env_ode$x')
        assert symbol_a.initial_value == 2.0
        assert symbol_t.initial_value is None
        assert symbol_x.initial_value is None
        assert symbol_y.initial_value is None
        assert symbol_x_orig.initial_value is None
        assert symbol_a.units == literals_model.units.get_unit('mV')
        assert symbol_t.units == literals_model.units.get_unit('ms')
        assert symbol_x.units == literals_model.units.get_unit('nA')
        assert symbol_y.units == literals_model.units.get_unit('per_pA')
        assert symbol_x_orig.units == literals_model.units.get_unit('pA')
        assert len(literals_model.equations) == 4
        assert str(literals_model.equations[0]) == 'Eq(Derivative(_env_ode$sv1, _environment$time), _1.0)'
        assert str(literals_model.equations[1]) == 'Eq(_env_ode$y, _1.0/_env_ode$x)'
        assert str(literals_model.equations[2]) == 'Eq(_env_ode$x_converted, 0.001*_1.0)'
        assert str(literals_model.equations[3]) == 'Eq(_env_ode$x, 1000.0*_env_ode$x_converted)'

    def test_add_input_free_variable_multiple(self, multiode_freevar_model):
        """ Tests the Model.convert_variable function that changes units of given variable.
        This particular case tests changing a free variable where there are multiple ode instances.
        For example::

            Original model
                var{time} time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};
                var y: mV {init: 3};

                ode(sv1, time) = 1{mV_per_ms};
                ode(y, time) = 2{mV_per_ms}

        convert_variable(time, second, DataDirectionFlow.INPUT)

            becomes
                var time: ms;
                var{time} time_converted: s;
                var{sv11} sv1: mV {init: 2};
                var sv1_orig_deriv mV_per_ms
                var y : mV {init: 3}
                var y_orig_deriv mV_per_ms

                time = 1000 * time_converted;
                sv1_orig_deriv = 1{mV_per_ms}
                ode(sv1, time_converted) = 1000 * sv1_orig_deriv
                y_orig_deriv = 2{mV_per_ms}
                ode(y, time_converted) = 1000 * y_orig_deriv
        """
        ms_unit = multiode_freevar_model.units.get_unit('ms')
        second_unit = multiode_freevar_model.units.get_unit('second')
        original_var = multiode_freevar_model.get_variable_by_name('environment$time')

        assert self._original_state_multiode_freevar(multiode_freevar_model)
        # test no change in units
        multiode_freevar_model.convert_variable(original_var, ms_unit, DataDirectionFlow.INPUT)
        assert self._original_state_multiode_freevar(multiode_freevar_model)

        # change ms to s
        multiode_freevar_model.convert_variable(original_var, second_unit, DataDirectionFlow.INPUT)
        assert len(multiode_freevar_model.variables()) == 7
        symbol_a = multiode_freevar_model.get_variable_by_cmeta_id('sv11')
        assert symbol_a.initial_value == 2.0
        assert symbol_a.units == multiode_freevar_model.units.get_unit('mV')
        assert symbol_a.name == 'env_ode$sv1'
        symbol_t = multiode_freevar_model.get_variable_by_cmeta_id('time')
        assert symbol_t.units == 'second'
        assert symbol_t.name == 'environment$time_converted'
        symbol_orig = multiode_freevar_model.get_variable_by_name('env_ode$sv1')
        mV = multiode_freevar_model.units.get_unit('mV')
        ms = multiode_freevar_model.units.get_unit('ms')
        assert symbol_orig.units == mV
        symbol_derv = multiode_freevar_model.get_variable_by_name('env_ode$sv1_orig_deriv')
        assert symbol_derv.units == mV / ms
        symbol_orig_y = multiode_freevar_model.get_variable_by_name('env_ode$y')
        assert symbol_orig_y.units == mV
        assert symbol_orig_y.initial_value == 3.0
        symbol_derv_y = multiode_freevar_model.get_variable_by_name('env_ode$y_orig_deriv')
        assert symbol_derv_y.units == mV / ms
        assert len(multiode_freevar_model.equations) == 5
        assert str(multiode_freevar_model.equations[0]) == 'Eq(_environment$time, 1000.0*_environment$time_converted)'
        assert str(multiode_freevar_model.equations[1]) == 'Eq(_env_ode$sv1_orig_deriv, _1.0)'
        assert str(multiode_freevar_model.equations[2]) == 'Eq(Derivative(_env_ode$sv1, ' \
                                                           '_environment$time_converted), ' \
                                                           '1000.0*_env_ode$sv1_orig_deriv)'
        assert str(multiode_freevar_model.equations[3]) == 'Eq(_env_ode$y_orig_deriv, _2.0)'
        assert str(multiode_freevar_model.equations[4]) == 'Eq(Derivative(_env_ode$y, _environment$time_converted), ' \
                                                           '1000.0*_env_ode$y_orig_deriv)'

    def test_multiple_odes(self, multiode_model):
        """ Tests the Model.convert_variable function that changes units of given variable.
        This tests that any other uses of the derivative on the rhs of other equations
        are replaced with the new variable representing the old derivative
        For example::

            Original model
                var time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};
                var x: mV_per_ms;
                var y: mV {init: 3};

                ode(sv1, time) = 1{mV_per_ms};
                x = 1{mV_per_ms} + ode(sv1, time) * 3{mV_per_ms};
                ode(y, time) = 2{mV_per_ms};

        convert_variable(sv1, volt, DataDirectionFlow.INPUT)

        becomes
                var{time} time: ms {pub: in};
                var{sv11} sv1_converted: V {init: 0.002};
                var sv1 mV {init: 2}
                var sv1_orig_deriv mV_per_ms
                var x: mV_per_ms;
                var y: mV {init: 3};
                var y_orig_deriv

                ode(y, time) = 2{mv_per_ms};
                sv1 = 1000 * sv1_converted
                sv1_orig_deriv = 1{mV_per_ms}
                ode(sv1_converted, time) = 0.001 * sv1_orig_deriv
                x = 1{mV_per_ms} + 3 * sv1_orig_deriv
        """
        mV_unit = multiode_model.units.get_unit('mV')
        volt_unit = multiode_model.units.get_unit('volt')
        original_var = multiode_model.get_variable_by_name('env_ode$sv1')

        assert self._original_state_multiode_model(multiode_model)
        # test no change in units
        multiode_model.convert_variable(original_var, mV_unit, DataDirectionFlow.INPUT)
        assert self._original_state_multiode_model(multiode_model)

        # change mV to V
        multiode_model.convert_variable(original_var, volt_unit, DataDirectionFlow.INPUT)
        assert len(multiode_model.equations) == 5
        assert str(multiode_model.equations[0]) == 'Eq(Derivative(_env_ode$y, _environment$time), _2.0)'
        assert str(multiode_model.equations[1]) == 'Eq(_env_ode$sv1, 1000.0*_env_ode$sv1_converted)'
        assert str(multiode_model.equations[2]) == 'Eq(_env_ode$sv1_orig_deriv, _1.0)'
        assert str(multiode_model.equations[3]) == 'Eq(Derivative(_env_ode$sv1_converted, _environment$time), ' \
                                                   '0.001*_env_ode$sv1_orig_deriv)'
        assert str(multiode_model.equations[4]) == 'Eq(_env_ode$x, _1.0 + _3.0*_env_ode$sv1_orig_deriv)'

    def test_multiple_odes_1(self, multiode_model):
        """ Tests the Model.convert_variable function that changes units of given variable.
        This tests that any other uses of the derivative on the rhs of other equations
        are replaced with the new variable representing the old derivative
        For example::

            Original model
                var time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};
                var x: mV_per_ms;
                var y: mV {init: 3};

                ode(sv1, time) = 1{mV_per_ms};
                x = 1{mV_per_ms} + ode(sv1, time) * 3{mV_per_ms};
                ode(y, time) = 2{mV_per_ms};


        convert_variable(time, second, DataDirectionFlow.INPUT)

            becomes
                var time: ms;
                var{time} time_converted: s;
                var{sv11} sv1: mV {init: 2};
                var sv1_orig_deriv mV_per_ms
                var x: mV_per_ms;
                var y: mV {init: 3};
                var y_orig_deriv: {mv_per_ms}

                time_converted = 1000 * time
                y_orig_deriv = 2{mV_per_ms}
                ode(y, time_converted) = 1000 * y_orig_deriv;
                sv1_orig_deriv = 1{mV_per_ms}
                ode(sv1, time_converted) = 1000 * sv1_orig_deriv
                x = 1 + 3 * sv1_orig_deriv
        """
        ms_unit = multiode_model.units.get_unit('ms')
        second_unit = multiode_model.units.get_unit('second')
        original_var = multiode_model.get_variable_by_name('environment$time')

        old_states = multiode_model.get_state_variables()
        assert len(old_states) == 2

        assert self._original_state_multiode_model(multiode_model)
        # test no change in units
        multiode_model.convert_variable(original_var, ms_unit, DataDirectionFlow.INPUT)
        assert self._original_state_multiode_model(multiode_model)

        # change ms to s
        multiode_model.convert_variable(original_var, second_unit, DataDirectionFlow.INPUT)
        assert len(multiode_model.equations) == 6
        assert str(multiode_model.equations[0]) == 'Eq(_environment$time, 1000.0*_environment$time_converted)'
        assert str(multiode_model.equations[1]) == 'Eq(_env_ode$sv1_orig_deriv, _1.0)'
        assert str(multiode_model.equations[2]) == 'Eq(Derivative(_env_ode$sv1, _environment$time_converted), ' \
                                                   '1000.0*_env_ode$sv1_orig_deriv)'
        assert str(multiode_model.equations[3]) == 'Eq(_env_ode$y_orig_deriv, _2.0)'
        assert str(multiode_model.equations[4]) == 'Eq(Derivative(_env_ode$y, _environment$time_converted), ' \
                                                   '1000.0*_env_ode$y_orig_deriv)'
        assert str(multiode_model.equations[5]) == 'Eq(_env_ode$x, _1.0 + _3.0*_env_ode$sv1_orig_deriv)'

        # test the graph forming and get_state_variables still works
        states = multiode_model.get_state_variables()
        assert len(states) == 2
        assert old_states == states

    def test_add_output(self, literals_model):
        """ Tests the Model.convert_variable function that changes units of given variable.
         This tests the case when variable to be changed is an OUTPUT

        For example::

            Original model
                var{time} time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};
                var{current} x: pA;
                var y: per_pA;

                ode(sv1, time) = 1{mV_per_ms};
                x = 1{pA};
                y = 1{dimensionless}/x;

        convert_variable(x, nA, DataDirectionFlow.OUTPUT)

        becomes
                var time: ms;
                var{sv11} sv1: mV {init: 2};
                var{current} x_converted: nA;
                var x : pA
                var y: per_pA;

                ode(sv1, time) = 1{mV_per_ms};
                x = 1{pA}
                y = 1{dimensionless}/x;
                x_converted = 0.001 * x
        """
        pA_unit = literals_model.units.get_unit('pA')
        nA_unit = literals_model.units.get_unit('nA')
        original_var = literals_model.get_variable_by_name('env_ode$x')

        assert self._original_state_literals_model(literals_model)
        # test no change in units
        literals_model.convert_variable(original_var, pA_unit, DataDirectionFlow.OUTPUT)
        assert self._original_state_literals_model(literals_model)

        # change pA to nA
        literals_model.convert_variable(original_var, nA_unit, DataDirectionFlow.OUTPUT)
        assert len(literals_model.variables()) == 6
        symbol_a = literals_model.get_variable_by_cmeta_id('sv11')
        symbol_t = literals_model.get_variable_by_cmeta_id('time')
        symbol_x = literals_model.get_variable_by_cmeta_id('current')
        assert symbol_x.name == 'env_ode$x_converted'
        symbol_y = literals_model.get_variable_by_name('env_ode$y')
        symbol_x_orig = literals_model.get_variable_by_name('env_ode$x')
        assert symbol_a.initial_value == 2.0
        assert symbol_t.initial_value is None
        assert symbol_x.initial_value is None
        assert symbol_y.initial_value is None
        assert symbol_x_orig.initial_value is None
        assert symbol_a.units == literals_model.units.get_unit('mV')
        assert symbol_t.units == literals_model.units.get_unit('ms')
        assert symbol_x.units == literals_model.units.get_unit('nA')
        assert symbol_y.units == literals_model.units.get_unit('per_pA')
        assert symbol_x_orig.units == literals_model.units.get_unit('pA')
        assert len(literals_model.equations) == 4
        assert str(literals_model.equations[0]) == 'Eq(Derivative(_env_ode$sv1, _environment$time), _1.0)'
        assert str(literals_model.equations[1]) == 'Eq(_env_ode$x, _1.0)'
        assert str(literals_model.equations[2]) == 'Eq(_env_ode$y, _1.0/_env_ode$x)'
        assert str(literals_model.equations[3]) == 'Eq(_env_ode$x_converted, 0.001*_env_ode$x)'
        state_variables = literals_model.get_state_variables()
        assert len(state_variables) == 1
        assert symbol_a in state_variables

    def test_add_output_state_variable(self, local_model):
        """ Tests the Model.convert_variable function that changes units of given variable.
        This particular test is when a state variable is being converted as an output

        For example::

            Original model
                var{time} time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};

                ode(sv1, time) = 1{mV_per_ms};

        convert sv11 from mV to V

            creates model
                var{time} time: ms {pub: in};
                var{sv11} sv1_converted: V;
                var sv1 mV {init: 2}

                ode(sv1, time) = 1{mV_per_ms};
                sv1_converted = sv1 / 1000
        """
        mV_unit = local_model.units.get_unit('mV')
        volt_unit = local_model.units.get_unit('volt')
        original_var = local_model.get_variable_by_name('env_ode$sv1')

        assert self._original_state_local_model(local_model)
        # test no change in units
        local_model.convert_variable(original_var, mV_unit, DataDirectionFlow.OUTPUT)
        assert self._original_state_local_model(local_model)

        # change mV to V
        local_model.convert_variable(original_var, volt_unit, DataDirectionFlow.OUTPUT)
        assert len(local_model.variables()) == 4
        symbol_a = local_model.get_variable_by_cmeta_id('sv11')
        assert symbol_a.initial_value is None
        assert symbol_a.units == local_model.units.get_unit('volt')
        assert symbol_a.name == 'env_ode$sv1_converted'
        symbol_t = local_model.get_variable_by_cmeta_id('time')
        assert symbol_t.units == local_model.units.get_unit('ms')
        symbol_orig = local_model.get_variable_by_name('env_ode$sv1')
        assert symbol_orig.units == local_model.units.get_unit('mV')
        assert symbol_orig.initial_value == 2.0
        assert len(local_model.equations) == 2
        assert str(local_model.equations[0]) == 'Eq(Derivative(_env_ode$sv1, _environment$time), _1.0)'
        assert str(local_model.equations[1]) == 'Eq(_env_ode$sv1_converted, 0.001*_env_ode$sv1)'
        state_variables = local_model.get_state_variables()
        assert len(state_variables) == 1
        assert symbol_orig in state_variables

    def test_add_output_free_variable(self, local_model):
        """ Tests the Model.convert_variable function that changes units of given variable.
        This particular test is when a free variable is being converted as an output

        For example::

            Original model
                var{time} time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};

                ode(sv1, time) = 1{mV_per_ms};

        convert time from ms to s

            creates model
                var{time} time: ms {pub: in};
                var{sv11} sv1: mV {init: 2};
                var{time} time_converted: s

                ode(sv1, time) = 1{mV_per_ms};
                time_converted = 0.001 * time
        """

        ms_unit = local_model.units.get_unit('ms')
        second_unit = local_model.units.get_unit('second')
        original_var = local_model.get_variable_by_name('environment$time')

        assert self._original_state_local_model(local_model)
        # test no change in units
        local_model.convert_variable(original_var, ms_unit, DataDirectionFlow.OUTPUT)
        assert self._original_state_local_model(local_model)

        # change ms to s
        local_model.convert_variable(original_var, second_unit, DataDirectionFlow.OUTPUT)
        assert len(local_model.variables()) == 4
        symbol_a = local_model.get_variable_by_cmeta_id('sv11')
        assert symbol_a.initial_value == 2.0
        assert symbol_a.units == local_model.units.get_unit('mV')
        assert symbol_a.name == 'env_ode$sv1'
        symbol_t = local_model.get_variable_by_cmeta_id('time')
        assert symbol_t.units == local_model.units.get_unit('second')
        assert symbol_t.name == 'environment$time_converted'
        symbol_orig = local_model.get_variable_by_name('environment$time')
        assert symbol_orig.units == local_model.units.get_unit('ms')
        assert len(local_model.equations) == 2
        assert str(local_model.equations[0]) == 'Eq(Derivative(_env_ode$sv1, _environment$time), _1.0)'
        assert str(local_model.equations[1]) == 'Eq(_environment$time_converted, 0.001*_environment$time)'
        state_variables = local_model.get_state_variables()
        assert len(state_variables) == 1
        assert symbol_a in state_variables
        assert local_model.get_free_variable() == symbol_orig

    def test_convert_variable_invalid_arguments(self, local_model):
        """ Tests the Model.convert_variable() function when involid arguments are passed.
        """
        unit = local_model.units.get_unit('second')
        variable = local_model.get_free_variable()
        direction = DataDirectionFlow.INPUT
        bad_unit = local_model.units.get_unit('mV')

        # arguments wrong types
        with pytest.raises(AssertionError):
            local_model.convert_variable('x', unit, direction)

        with pytest.raises(AssertionError):
            local_model.convert_variable(variable, 'x', direction)

        with pytest.raises(AssertionError):
            local_model.convert_variable(variable, unit, 'x')

        # variable not present in model
        model = shared.load_model('literals_for_conversion_tests')
        other_var = model.get_variable_by_name('env_ode$x')
        with pytest.raises(AssertionError):
            local_model.convert_variable(other_var, unit, direction)

        # ontology term not present in model
        with pytest.raises(AssertionError):
            local_model.convert_variable('current', unit, direction)

        # unit conversion is impossible
        with pytest.raises(DimensionalityError):
            local_model.convert_variable(variable, bad_unit, direction)

    def test_convert_same_unit_different_name(self, br_model):
        """ Tests the Model.convert_variable() function when conversion to current unit under a different name."""
        br_model.units.add_unit('millimolar', 'mole / 1000 / litre')
        unit = br_model.units.get_unit('concentration_units')
        variable = br_model.get_variable_by_ontology_term((shared.OXMETA, "cytosolic_calcium_concentration"))
        direction = DataDirectionFlow.INPUT
        assert br_model.convert_variable(variable, unit, direction) == variable

    def test_convert_variable_no_cmeta_id(self, br_model):
        """
        Tests converting units works on a variable without a cmeta id (there was a bug before that stopped this from
        working).
        """
        time = br_model.get_free_variable()
        time_units = br_model.units.evaluate_units(time)
        assert br_model.units.format(time_units) == 'ms'
        desired_units = br_model.units.get_unit('second')
        converted_time = br_model.convert_variable(time, desired_units, DataDirectionFlow.INPUT)
        converted_time_units = br_model.units.evaluate_units(converted_time)
        assert br_model.units.format(converted_time_units) == 'second'

    def test_unique_names(self, silly_names):
        # original state
        def test_original_state(silly_names):
            assert len(silly_names.variables()) == 5
            symbol_a = silly_names.get_variable_by_cmeta_id('sv11')
            symbol_t = silly_names.get_variable_by_cmeta_id('time')
            assert symbol_a.initial_value == 2.0
            assert symbol_a.units == silly_names.units.get_unit('mV')
            assert symbol_t.units == silly_names.units.get_unit('ms')
            assert silly_names.get_variable_by_name('env_ode$sv1_converted')
            assert silly_names.get_variable_by_name('env_ode$sv1_orig_deriv')
            assert len(silly_names.equations) == 1
            assert str(silly_names.equations[0]) == 'Eq(Derivative(_env_ode$sv1, _environment$time), _1.0)'
            state_variables = silly_names.get_state_variables()
            assert len(state_variables) == 1
            assert symbol_a in state_variables
            assert silly_names.get_free_variable() == symbol_t
            return True

        volt_unit = silly_names.units.get_unit('volt')
        original_var = silly_names.get_variable_by_name('env_ode$sv1')

        assert test_original_state(silly_names)
        # change mV to V
        silly_names.convert_variable(original_var, volt_unit, DataDirectionFlow.INPUT)
        assert len(silly_names.variables()) == 7
        symbol_a = silly_names.get_variable_by_cmeta_id('sv11')
        assert symbol_a.initial_value == 0.002
        assert symbol_a.units == silly_names.units.get_unit('volt')
        assert symbol_a.name == 'env_ode$sv1_converted_a'
        symbol_t = silly_names.get_variable_by_cmeta_id('time')
        assert symbol_t.units == silly_names.units.get_unit('ms')
        assert symbol_t.name == 'environment$time'
        symbol_orig = silly_names.get_variable_by_name('env_ode$sv1')
        assert symbol_orig.units == silly_names.units.get_unit('mV')
        assert silly_names.get_variable_by_name('env_ode$sv1_converted')
        assert silly_names.get_variable_by_name('env_ode$sv1_orig_deriv')
        symbol_derv = silly_names.get_variable_by_name('env_ode$sv1_orig_deriv_a')
        assert symbol_derv.units == silly_names.units.get_unit('mV') / silly_names.units.get_unit('ms')
        assert symbol_derv.initial_value is None
        assert len(silly_names.equations) == 3
        assert str(silly_names.equations[0]) == 'Eq(_env_ode$sv1, 1000.0*_env_ode$sv1_converted_a)'
        assert str(silly_names.equations[1]) == 'Eq(_env_ode$sv1_orig_deriv_a, _1.0)'
        assert str(silly_names.equations[2]) == 'Eq(Derivative(_env_ode$sv1_converted_a, _environment$time), ' \
                                                '0.001*_env_ode$sv1_orig_deriv_a)'
