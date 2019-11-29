from collections import OrderedDict

import pytest
import sympy as sp

from cellmlmanip.model import NumberDummy, VariableDummy
from cellmlmanip.units import (
    BooleanUnitsError,
    ExpressionWithUnitPrinter,
    InputArgumentMustBeNumberError,
    InputArgumentsInvalidUnitsError,
    InputArgumentsMustBeDimensionlessError,
    UnexpectedMathUnitsError,
    UnitCalculator,
    UnitStore,
)


class TestUnits(object):
    # These represent CellML <units><unit>...</unit></units> elements
    test_definitions = OrderedDict({
        'ms': [{'units': 'second', 'prefix': 'milli'}],
        'usec': [{'units': 'second', 'prefix': 'micro'}],
        'mV': [{'units': 'volt', 'prefix': 'milli'}],
        'uV': [{'units': 'volt', 'prefix': 'micro'}],
        'mM': [{'prefix': 'milli', 'units': 'mole'}, {'units': 'litre', 'exponent': '-1'}],
        'milli_mole': [{'prefix': 'milli', 'units': 'mole'}],
        'millisecond': [{'prefix': 'milli', 'units': 'second'}],
    })
    test_definitions.update({
        'per_ms': [{'units': 'ms', 'exponent': '-1'}],
        'per_mV': [{'units': 'volt', 'prefix': 'milli', 'exponent': '-1'}],
        'mV_per_ms': [{'units': 'mV', 'exponent': '1'}, {'units': 'ms', 'exponent': '-1'}],
        'mV_per_s': [{'units': 'mV', 'exponent': '1'}, {'units': 'second', 'exponent': '-1'}],
        'mV_per_usec': [
            {'units': 'mV', 'exponent': '1'}, {'prefix': 'micro', 'units': 'second', 'exponent': '-1'}],
        'mM_per_ms': [{'units': 'mM'}, {'units': 'ms', 'exponent': '-1'}],
        'ms_power_prefix': [{'prefix': '-3', 'units': 'second'}],
        'ms_with_multiplier': [{'multiplier': 0.001, 'units': 'second'}],
    })

    @pytest.fixture(scope="class")
    def quantity_store(self):
        qs = UnitStore(model=None)
        for unit_name, unit_attributes in self.test_definitions.items():
            qs.add_custom_unit(unit_name, unit_attributes)
        return qs

    def test_quantity_translation(self, quantity_store):
        unit_registry = quantity_store.ureg

        assert quantity_store.get_quantity('dimensionless') == unit_registry.dimensionless

        # Units defined in the test CellML <model>:
        for name in TestUnits.test_definitions.keys():
            assert quantity_store.get_quantity(name) != unit_registry.Quantity(1)

        # Pint built-in units
        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('kilogram'),
            unit_registry.kilogram
        )

        # Custom units defined in CellML example
        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('per_ms'),
            unit_registry.millisecond**(-1)
        )

        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('mM_per_ms'),
            (unit_registry.milli_mole / unit_registry.liter) / unit_registry.ms
        )

        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('mV_per_usec'),
            unit_registry.mV / unit_registry.usec
        )

        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('ms_power_prefix'),
            unit_registry.millisecond
        )

        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('ms_with_multiplier'),
            unit_registry.millisecond
        )

    def test_conversion_factor(self, quantity_store):
        ureg = quantity_store.ureg
        assert quantity_store.get_conversion_factor(quantity=1 * ureg.ms, to_unit=ureg.second) == 0.001
        assert quantity_store.get_conversion_factor(quantity=1 * ureg.volt, to_unit=ureg.mV) == 1000.0

        assert quantity_store.get_conversion_factor(
            quantity=1 * quantity_store.get_quantity('milli_mole'),
            to_unit=quantity_store.get_quantity('mole')
        ) == 0.001

    def test_add_custom_unit(self):
        unitstore = UnitStore(model=None)
        assert (unitstore._is_unit_defined('newunit') is False)
        # add a new unit called 'newunit' which consists of
        # (2 * second) - I know not realistic
        unit_attributes = [{'multiplier': 2, 'units': 'second'}]
        unitstore.add_custom_unit('newunit', unit_attributes)
        assert (unitstore._is_unit_defined('newunit') is True)

    def test_add_custom_unit_1(self):
        unitstore = UnitStore(model=None)
        assert (unitstore._is_unit_defined('newunit') is False)
        # add a new unit called 'newunit' which consists of
        # (millimetres)
        unit_attributes = [{'prefix': 'milli', 'units': 'metre'}]
        unitstore.add_custom_unit('newunit', unit_attributes)
        assert (unitstore._is_unit_defined('newunit') is True)

    def test_add_custom_unit_2(self):
        unitstore = UnitStore(model=None)
        assert (unitstore._is_unit_defined('newunit') is False)
        # add a new unit called 'newunit' which consists of
        # (milliVolts)
        unit_attributes = [{'prefix': -3, 'units': 'volt'}]
        unitstore.add_custom_unit('newunit', unit_attributes)
        assert (unitstore._is_unit_defined('newunit') is True)

    def test_add_custom_unit_3(self):
        unitstore = UnitStore(model=None)
        assert (unitstore._is_unit_defined('newunit') is False)
        # add a new unit called 'newunit' which consists of
        # (metre per second)
        unit_attributes = [{'units': 'metre'}, {'units': 'second', 'exponent': -1}]
        unitstore.add_custom_unit('newunit', unit_attributes)
        assert (unitstore._is_unit_defined('newunit') is True)

    def test_add_custom_unit_existing(self):
        unitstore = UnitStore(model=None)
        unit_attributes = [{'multiplier': 1, 'units': 'second'}]
        with pytest.raises(AssertionError):
            unitstore.add_custom_unit('second', unit_attributes)
            pytest.fail("Cannot redefine CellML unit <second>")
            pass

    def test_is_unit_defined(self, quantity_store):
        assert (quantity_store._is_unit_defined('ms') is True)
        assert (quantity_store._is_unit_defined('not_a_unit') is False)

    def test_make_pint_unit_definition(self, quantity_store):
        unit_attributes = [{'prefix': -3, 'units': 'metre'},
                           {'prefix': 'milli', 'exponent': -1, 'units': 'second'},
                           {'multiplier': 2, 'units': 'kilograms'}
                           ]
        assert(quantity_store._make_pint_unit_definition('kg_mm_per_ms', unit_attributes) ==
               'kg_mm_per_ms=(meter * 1e-3)*(((second * 0.001))**-1)*(2 * kilogram)')

    # Test UnitCalculator class

    def test_unit_calculator(self, quantity_store):
        """ Tests the units.UnitCalculator class. """

        ureg = quantity_store.ureg

        a = VariableDummy('a', ureg.meter)
        b = VariableDummy('b', ureg.second)
        c = VariableDummy('c', ureg.gram)
        d = VariableDummy('d', ureg.meter)
        x = VariableDummy('x', ureg.kilogram)
        y = VariableDummy('y', ureg.volt)
        n = VariableDummy('n', ureg.dimensionless)
        _1 = NumberDummy(1, ureg.kelvin)
        _2 = NumberDummy(2, ureg.dimensionless)
        _25 = NumberDummy(2.5, ureg.dimensionless)
        av = VariableDummy('av', ureg.meter, initial_value=2)

        unit_calculator = UnitCalculator(ureg)

        # units of numbers
        assert unit_calculator.traverse(_1).units == ureg.kelvin
        assert unit_calculator.traverse(sp.sympify("1.0")).units == ureg.dimensionless
        assert unit_calculator.traverse(_2).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.sympify("2")).units == ureg.dimensionless
        assert unit_calculator.traverse(_25).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.sympify("2.5")).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.sympify("12")).units == ureg.dimensionless

        # functions were all args should have same unit
        # and this is the unit that will be returned
        # pass cases
        assert unit_calculator.traverse(-a).units == ureg.meter
        assert unit_calculator.traverse(a + a + a + a).units == ureg.meter
        assert unit_calculator.traverse(a + 2*a + 3*a + 4*a).units == ureg.meter  # noqa: E226
        assert unit_calculator.traverse(2 * a - a).units == ureg.meter
        assert unit_calculator.traverse(a + d).units == ureg.meter
        assert unit_calculator.traverse(a - d).units == ureg.meter
        assert unit_calculator.traverse(sp.Abs(-2 * y)).units == ureg.volt
        assert unit_calculator.traverse(sp.floor(x)).units == ureg.kilogram
        result = unit_calculator.traverse(sp.floor(12.5) * a)
        assert result.units == ureg.meter and result.magnitude == 12.0 * a
        assert unit_calculator.traverse(sp.floor(_1)).units == ureg.kelvin
        result = unit_calculator.traverse(sp.floor(_25))
        assert result.units == ureg.dimensionless and result.magnitude == 2.0
        assert unit_calculator.traverse(sp.ceiling(x)).units == ureg.kilogram
        result = unit_calculator.traverse(sp.ceiling(12.6) * a)
        assert result.units == ureg.meter and result.magnitude == 13.0 * a
        assert unit_calculator.traverse(sp.ceiling(_1)).units == ureg.kelvin
        result = unit_calculator.traverse(sp.ceiling(_25))
        assert result.units == ureg.dimensionless and result.magnitude == 3.0

        # functions were all args should have same unit
        # and this is the unit that will be returned
        # fail cases
        with pytest.raises(InputArgumentsInvalidUnitsError):
            unit_calculator.traverse(a + b)
        try:
            unit_calculator.traverse(a + b)
        except InputArgumentsInvalidUnitsError as err:
            assert err.message == 'The arguments to this expression should all have the same units.'
            assert err.expression == '_a + _b'
        with pytest.raises(InputArgumentsInvalidUnitsError):
            unit_calculator.traverse(a - b)
        try:
            unit_calculator.traverse(a - b)
        except InputArgumentsInvalidUnitsError as err:
            assert err.message == 'The arguments to this expression should all have the same units.'
            assert err.expression == '_a - _b'
        try:
            sp.floor(x, y)
        except TypeError as err:
            assert err.args[0] == 'floor takes exactly 1 argument (2 given)'

        # special case - piecewise
        assert unit_calculator.traverse(sp.Piecewise((a, x < 1),
                                                     (a + a, x > 1),
                                                     (3 * a, True))).units == ureg.meter
        # fail special case -piecewise
        with pytest.raises(InputArgumentsInvalidUnitsError):
            unit_calculator.traverse(sp.Piecewise((a, x < 1), (b, x > 1), (c, True)))
        try:
            unit_calculator.traverse(sp.Piecewise((a, x < 1), (b, x > 1), (c, True)))
        except InputArgumentsInvalidUnitsError as err:
            assert err.message == 'The arguments to this expression should all have the same units.'
            assert err.expression == 'Piecewise((_a, _x < 1), (_b, _x > 1), (_c, True))'

        # cases with any units allowed as arguments
        assert unit_calculator.traverse((a * a) / b).units == ureg.meter**2 / ureg.second

        # root and power
        assert unit_calculator.traverse(a**_2).units == ureg.meter**2
        assert unit_calculator.traverse(sp.sqrt(c ** 2)).units == ureg.gram
        assert unit_calculator.traverse(sp.sqrt(a * d)).units == ureg.meter
        # root and power fails
        expr = a ** _1
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            unit_calculator.traverse(expr)
        try:
            unit_calculator.traverse(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The second argument to this expression should be dimensionless.'
            assert err.expression == '_a**_1'
        expr = a + a ** n
        with pytest.raises(InputArgumentMustBeNumberError):
            unit_calculator.traverse(expr)
        try:
            unit_calculator.traverse(expr)
        except InputArgumentMustBeNumberError as err:
            assert err.message == 'The second argument to this expression should be a number.'
            assert err.expression == '_a**_n'

        # bizarre cases
        expr = av ** _25
        result = unit_calculator.traverse(expr)
        assert result.units == ureg.meter**2.5  # and result.magnitude == 5.6568
        expr = sp.root(av, _25)
        result = unit_calculator.traverse(expr)
        assert result.units == ureg.meter**0.4  # and result.magnitude == 1.319507

        # relational operators throw an exception
        expr = a > _1
        with pytest.raises(BooleanUnitsError):
            unit_calculator.traverse(expr)

        # special case - derivative
        dadb = sp.diff(a * b, b)
        assert unit_calculator.traverse(dadb).units == ureg.meter

        dadb = sp.Derivative(a * b, b)
        assert unit_calculator.traverse(dadb).units == ureg.meter

        # log and exponent
        assert unit_calculator.traverse(sp.log(_2)).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.log(n)).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.ln(_2)).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.log(n, _2)).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.log(_2, 3)).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.exp(_2)).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.exp(n)).units == ureg.dimensionless
        # log and exponent - fails
        expr = sp.log(a)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            unit_calculator.traverse(expr)
        try:
            unit_calculator.traverse(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The arguments to this expression should be dimensionless.'
            assert err.expression == 'log(_a)'
        expr = sp.log(_2, av)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            unit_calculator.traverse(expr)
        try:
            unit_calculator.traverse(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The arguments to this expression should be dimensionless.'
            assert err.expression == 'log(_av)'
        expr = sp.exp(av)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            unit_calculator.traverse(expr)
        try:
            unit_calculator.traverse(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The arguments to this expression should be dimensionless.'
            assert err.expression == 'exp(_av)'

        # trig functions
        assert unit_calculator.traverse(sp.sin(n)).units == ureg.dimensionless
        # trig functions - fails
        expr = sp.cos(av)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            unit_calculator.traverse(expr)
        try:
            unit_calculator.traverse(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The arguments to this expression should be dimensionless.'
            assert err.expression == 'cos(_av)'

        # constants
        assert unit_calculator.traverse(sp.pi).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.E).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.oo).units == ureg.dimensionless
        assert unit_calculator.traverse(sp.nan).units == ureg.dimensionless

        expr = sp.true
        with pytest.raises(BooleanUnitsError):
            unit_calculator.traverse(expr)
        try:
            unit_calculator.traverse(expr)
        except BooleanUnitsError as err:
            assert err.message == 'This expression involves boolean values which do not conform to ' \
                                  'unit dimensionality rules.'
            assert err.expression == 'True'

        # factorial needs its own catch
        assert unit_calculator.traverse(sp.factorial(n)).units == ureg.dimensionless
        expr = sp.factorial(av)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            unit_calculator.traverse(expr)
        try:
            unit_calculator.traverse(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The arguments to this expression should be dimensionless.'
            assert err.expression == 'factorial(_av)'

        # logic functions throw an exception
        expr = a & b
        with pytest.raises(BooleanUnitsError):
            unit_calculator.traverse(expr)
        try:
            unit_calculator.traverse(expr)
        except BooleanUnitsError as err:
            assert err.message == 'This expression involves boolean values which do not conform to ' \
                                  'unit dimensionality rules.'
            assert err.expression == '_a & _b'

        # check that not gets caught as it is listed separately
        # does indeed get caught by is_Boolean
        expr = ~a
        with pytest.raises(BooleanUnitsError):
            unit_calculator.traverse(expr)

        # this is a generic test of a function with one dimensionless argument
        # this uses a function that is not in cellml - just to check the logic
        assert unit_calculator.traverse(sp.sign(n)).units == ureg.dimensionless

        expr = sp.Matrix([[1, 0], [0, 1]])
        with pytest.raises(UnexpectedMathUnitsError):
            unit_calculator.traverse(expr)

        N = sp.Matrix([[1, 0], [0, 1]])
        M = sp.Matrix([[1, 0], [0, 1]])
        expr = M + N
        with pytest.raises(UnexpectedMathUnitsError):
            unit_calculator.traverse(expr)

        expr = sp.cos(x).series(x, 0, 10)
        with pytest.raises(UnexpectedMathUnitsError):
            unit_calculator.traverse(expr)

    def test_expression_printer(self, quantity_store):
        """Tests the units.ExpressionWithUnitPrinter class. """

        ureg = quantity_store.ureg
        a = VariableDummy('a', ureg.meter)
        b = VariableDummy('b', ureg.second)
        y = VariableDummy('y', ureg.volt)
        _2 = NumberDummy(2, ureg.dimensionless)

        printer = ExpressionWithUnitPrinter()
        assert printer.doprint(a * a) == 'a[meter]**2'
        assert printer.doprint(a / b) == 'a[meter]/b[second]'
        assert printer.doprint(_2 * y) == '2.000000[dimensionless]*y[volt]'
        assert printer.doprint(sp.Derivative(a, b)) == 'Derivative(a[meter], b[second])'

    def test_dimensionally_equivalent(self, hh_model, OXMETA):
        membrane_stimulus_current_offset = hh_model.get_symbol_by_ontology_term(OXMETA,
                                                                                "membrane_stimulus_current_offset")
        membrane_stimulus_current_period = hh_model.get_symbol_by_ontology_term(OXMETA,
                                                                                "membrane_stimulus_current_period")
        membrane_voltage = hh_model.get_symbol_by_ontology_term(OXMETA, "membrane_voltage")
        assert hh_model.units.dimensionally_equivalent(membrane_stimulus_current_offset,
                                                       membrane_stimulus_current_period)
        assert not hh_model.units.dimensionally_equivalent(membrane_stimulus_current_offset, membrane_voltage)
