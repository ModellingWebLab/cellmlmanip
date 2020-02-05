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
from . import shared


class TestUnits(object):
    # These represent CellML <units><unit>...</unit></units> elements
    test_definitions = OrderedDict({
        'ms': 'second / 1000',
        'usec': 'second / 1000000',
        'mV': 'volt / 1000',
        'uV': 'volt / 1000000',
        'mM': '(mole / 1000) / litre',
        'milli_mole': 'mole / 1000',
        'millisecond': 'second / 1000',
    })
    test_definitions.update({
        'per_ms': '1 / ms',
        'per_mV': '1 / mV',
        'mV_per_ms': 'mV / ms',
        'mV_per_s': 'mV / second',
        'mV_per_usec': 'mV / usec',
        'mM_per_ms': 'mM / ms',
    })

    @pytest.fixture(scope="class")
    def unit_store(self):
        """A UnitStore for testing."""
        s = UnitStore()
        for name, definition in self.test_definitions.items():
            s.add_unit(name, definition)
        return s

    def test_quantity_translation(self, unit_store):
        """Tests that unit equality is correctly determined."""

        # Units defined in the test CellML <model>:
        for name in TestUnits.test_definitions.keys():
            assert unit_store.is_unit_defined(name)

        # Custom units defined in CellML example
        '''
        assert

        a = 1 * unit_store.get_unit('second')
        b = 1 * unit_store.get_unit('ms')

        print(a.to_base_units())
        print(b.to_base_units())




        assert unit_store.is_equal(
            1 * unit_store.get_unit('per_ms'),
            1 / (unit_store.get_unit('second') / 1000)
        )

        assert unit_store.is_equal(
            unit_store.get_unit('mM_per_ms'),
            (unit_registry.milli_mole / unit_registry.liter) / unit_registry.ms
        )

        assert unit_store.is_equal(
            unit_store.get_unit('mV_per_usec'),
            unit_registry.mV / unit_registry.usec
        )

        assert unit_store.is_equal(
            unit_store.get_unit('ms_power_prefix'),
            unit_registry.millisecond
        )

        assert unit_store.is_equal(
            unit_store.get_unit('ms_with_multiplier'),
            unit_registry.millisecond
        )
        '''

    def test_conversion_factor(self, unit_store):
        """ Tests Units.get_conversion_factor() function. """
        assert unit_store.get_conversion_factor(
            quantity=1 * unit_store.get_unit('ms'),
            to_unit=unit_store.get_unit('second')) == 0.001
        assert unit_store.get_conversion_factor(
            quantity=1 * unit_store.get_unit('volt'),
            to_unit=unit_store.get_unit('mV')) == 1000.0

        assert unit_store.get_conversion_factor(
            quantity=1 * unit_store.get_unit('milli_mole'),
            to_unit=unit_store.get_unit('mole')
        ) == 0.001

    def test_add_unit_0(self):
        """Tests Units.add_unit()."""
        unitstore = UnitStore()
        assert (unitstore.is_unit_defined('newunit') is False)

        # Add a new unit called 'newunit' which is 2 * second
        unitstore.add_unit('newunit', 'second * 2')
        assert (unitstore.is_unit_defined('newunit') is True)

    def test_add_unit_1(self):
        """Tests Units.add_unit()."""
        unitstore = UnitStore()
        assert (unitstore.is_unit_defined('newunit') is False)

        # Add a new unit called 'newunit' which is millimetres
        unitstore.add_unit('newunit', 'metre / 1000')
        assert (unitstore.is_unit_defined('newunit') is True)

    def test_add_unit_3(self):
        """ Tests Units.add_unit() function. """
        unitstore = UnitStore()
        assert (unitstore.is_unit_defined('newunit') is False)

        # Add a new unit called 'newunit' which is metre per second
        unitstore.add_unit('newunit', 'meter / second')
        assert (unitstore.is_unit_defined('newunit') is True)

    def test_add_unit_existing(self):
        """ Tests exception for Units.add_unit() function when unit already exists. """
        unitstore = UnitStore()
        with pytest.raises(AssertionError):
            unitstore.add_unit('second', 'second * 2')
            pytest.fail('Cannot redefine CellML unit <second>')

    # def test_is_unit_defined(self, unit_store):
    #     """ Tests Units._is_unit_defined() function. """
    #     assert (unit_store.is_unit_defined('ms') is True)
    #    assert (unit_store.is_unit_defined('not_a_unit') is False)

    # def test_make_pint_unit_definition(self, unit_store):
    #    """ Tests Units._make_pint_unit_definition() function. """
    #    unit_attributes = [{'prefix': -3, 'units': 'metre'},
    #                       {'prefix': 'milli', 'exponent': -1, 'units': 'second'},
    #                       {'multiplier': 2, 'units': 'kilograms'}
    #                       ]
    #    assert(unit_store._make_pint_unit_definition('kg_mm_per_ms', unit_attributes) ==
    #           'kg_mm_per_ms=(meter * 1e-3)*(((second * 0.001))**-1)*(2 * kilogram)')

    '''
    def test_shared_unit_registry(self):
        """ Tests sharing a single unit registry between two UnitStores. """

        # Without shared registry
        store1 = UnitStore()
        store2 = UnitStore()
        assert not store1.is_unit_defined('mmm')
        assert not store2.is_unit_defined('mmm')

        # ...adding unit to store 1 doesn't affect store 2
        unit_attributes = [{'prefix': 'milli', 'units': 'metre'}]
        store1.add_custom_unit('mmm', unit_attributes)
        assert store1.is_unit_defined('mmm')
        assert not store2.is_unit_defined('mmm')
        store2.add_custom_unit('mmm', unit_attributes)
        assert store2.is_unit_defined('mmm')

        # ...and the units do not count as equivalent, because they're from different registries (even thought they have
        # the same basi SI unit expansion).
        with pytest.raises(ValueError):
            assert store1.ureg.mmm != store2.ureg.mmm

        # ...but quantities are the same
        assert store1.ureg.parse_expression('mmm') == store2.ureg.parse_expression('mmm')

        # With a shared registry
        store2 = UnitStore(store1.ureg)
        assert store1.is_unit_defined('mmm')
        assert not store2.is_unit_defined('mmm')
    '''

    # Test UnitCalculator class

    def test_unit_calculator(self, unit_store):
        """ Tests the units.UnitCalculator class. """

        a = VariableDummy('a', unit_store.get_unit('meter'))
        b = VariableDummy('b', unit_store.get_unit('second'))
        c = VariableDummy('c', unit_store.get_unit('gram'))
        d = VariableDummy('d', unit_store.get_unit('meter'))
        x = VariableDummy('x', unit_store.get_unit('kilogram'))
        y = VariableDummy('y', unit_store.get_unit('volt'))
        n = VariableDummy('n', unit_store.get_unit('dimensionless'))
        _1 = NumberDummy(1, unit_store.get_unit('kelvin'))
        _2 = NumberDummy(2, unit_store.get_unit('dimensionless'))
        _25 = NumberDummy(2.5, unit_store.get_unit('dimensionless'))
        av = VariableDummy('av', unit_store.get_unit('meter'), initial_value=2)

        unit_calculator = UnitCalculator(unit_store)

        # units of numbers
        assert unit_calculator.traverse(_1).units == unit_store.get_unit('kelvin')
        assert unit_calculator.traverse(sp.sympify("1.0")).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(_2).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.sympify("2")).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(_25).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.sympify("2.5")).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.sympify("12")).units == unit_store.get_unit('dimensionless')

        # functions were all args should have same unit
        # and this is the unit that will be returned
        # pass cases
        assert unit_calculator.traverse(-a).units == unit_store.get_unit('meter')
        assert unit_calculator.traverse(a + a + a + a).units == unit_store.get_unit('meter')
        assert unit_calculator.traverse(a + 2 * a + 3 * a + 4 * a).units == unit_store.get_unit('meter')
        assert unit_calculator.traverse(2 * a - a).units == unit_store.get_unit('meter')
        assert unit_calculator.traverse(a + d).units == unit_store.get_unit('meter')
        assert unit_calculator.traverse(a - d).units == unit_store.get_unit('meter')
        assert unit_calculator.traverse(sp.Abs(-2 * y)).units == unit_store.get_unit('volt')
        assert unit_calculator.traverse(sp.floor(x)).units == unit_store.get_unit('kilogram')
        result = unit_calculator.traverse(sp.floor(12.5) * a)
        assert result.units == unit_store.get_unit('meter') and result.magnitude == 12 * a
        assert unit_calculator.traverse(sp.floor(_1)).units == unit_store.get_unit('kelvin')
        result = unit_calculator.traverse(sp.floor(_25))
        assert result.units == unit_store.get_unit('dimensionless') and result.magnitude == 2.0
        assert unit_calculator.traverse(sp.ceiling(x)).units == unit_store.get_unit('kilogram')
        result = unit_calculator.traverse(sp.ceiling(12.6) * a)
        assert result.units == unit_store.get_unit('meter') and result.magnitude == 13 * a
        assert unit_calculator.traverse(sp.ceiling(_1)).units == unit_store.get_unit('kelvin')
        result = unit_calculator.traverse(sp.ceiling(_25))
        assert result.units == unit_store.get_unit('dimensionless') and result.magnitude == 3.0

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
                                                     (3 * a, True))).units == unit_store.get_unit('meter')
        # fail special case -piecewise
        with pytest.raises(InputArgumentsInvalidUnitsError):
            unit_calculator.traverse(sp.Piecewise((a, x < 1), (b, x > 1), (c, True)))
        try:
            unit_calculator.traverse(sp.Piecewise((a, x < 1), (b, x > 1), (c, True)))
        except InputArgumentsInvalidUnitsError as err:
            assert err.message == 'The arguments to this expression should all have the same units.'
            assert err.expression == 'Piecewise((_a, _x < 1), (_b, _x > 1), (_c, True))'

        # cases with any units allowed as arguments
        assert unit_calculator.traverse((a * a) / b).units == (
            unit_store.get_unit('meter') ** 2 / unit_store.get_unit('second'))

        # root and power
        assert unit_calculator.traverse(a**_2).units == unit_store.get_unit('meter')**2
        assert unit_calculator.traverse(sp.sqrt(c ** 2)).units == unit_store.get_unit('gram')
        assert unit_calculator.traverse(sp.sqrt(a * d)).units == unit_store.get_unit('meter')
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
        assert result.units == unit_store.get_unit('meter')**2.5  # and result.magnitude == 5.6568
        expr = sp.root(av, _25)
        result = unit_calculator.traverse(expr)
        assert result.units == unit_store.get_unit('meter')**0.4  # and result.magnitude == 1.319507

        # relational operators throw an exception
        expr = a > _1
        with pytest.raises(BooleanUnitsError):
            unit_calculator.traverse(expr)

        # special case - derivative
        dadb = sp.diff(a * b, b)
        assert unit_calculator.traverse(dadb).units == unit_store.get_unit('meter')

        dadb = sp.Derivative(a * b, b)
        assert unit_calculator.traverse(dadb).units == unit_store.get_unit('meter')

        # log and exponent
        assert unit_calculator.traverse(sp.log(_2)).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.log(n)).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.ln(_2)).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.log(n, _2)).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.log(_2, 3)).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.exp(_2)).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.exp(n)).units == unit_store.get_unit('dimensionless')
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
        assert unit_calculator.traverse(sp.sin(n)).units == unit_store.get_unit('dimensionless')
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
        assert unit_calculator.traverse(sp.pi).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.E).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.oo).units == unit_store.get_unit('dimensionless')
        assert unit_calculator.traverse(sp.nan).units == unit_store.get_unit('dimensionless')

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
        assert unit_calculator.traverse(sp.factorial(n)).units == unit_store.get_unit('dimensionless')
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
        assert unit_calculator.traverse(sp.sign(n)).units == unit_store.get_unit('dimensionless')

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

    def test_expression_printer(self, unit_store):
        """Tests the units.ExpressionWithUnitPrinter class. """

        a = VariableDummy('a', 'meter')
        b = VariableDummy('b', 'second')
        y = VariableDummy('y', 'volt')
        _2 = NumberDummy(2, 'dimensionless')

        printer = ExpressionWithUnitPrinter()
        assert printer.doprint(a * a) == 'a[meter]**2'
        assert printer.doprint(a / b) == 'a[meter]/b[second]'
        assert printer.doprint(_2 * y) == '2.000000[dimensionless]*y[volt]'
        assert printer.doprint(sp.Derivative(a, b)) == 'Derivative(a[meter], b[second])'

    def test_dimensionally_equivalent(self, hh_model):
        """ Tests Units.dimensionally_equivalent() function. """
        membrane_stimulus_current_offset = hh_model.get_symbol_by_ontology_term(shared.OXMETA,
                                                                                "membrane_stimulus_current_offset")
        membrane_stimulus_current_period = hh_model.get_symbol_by_ontology_term(shared.OXMETA,
                                                                                "membrane_stimulus_current_period")
        membrane_voltage = hh_model.get_symbol_by_ontology_term(shared.OXMETA, "membrane_voltage")
        assert hh_model.units.dimensionally_equivalent(membrane_stimulus_current_offset,
                                                       membrane_stimulus_current_period)
        assert not hh_model.units.dimensionally_equivalent(membrane_stimulus_current_offset, membrane_voltage)
