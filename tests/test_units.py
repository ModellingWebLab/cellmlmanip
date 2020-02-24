import pint
import pytest
import sympy as sp

from cellmlmanip.model import NumberDummy, VariableDummy
from cellmlmanip.units import (
    BooleanUnitsError,
    InputArgumentMustBeNumberError,
    InputArgumentsInvalidUnitsError,
    InputArgumentsMustBeDimensionlessError,
    UnexpectedMathUnitsError,
    UnitStore,
)


class TestUnits(object):

    def test_add_unit(self):
        """Tests UnitStore.add_unit()."""

        # Add ordinary unit
        unitstore = UnitStore()
        assert not unitstore.is_defined('u1')
        u1 = unitstore.add_unit('u1', 'second * 2')
        assert unitstore.is_defined('u1')
        assert u1 == unitstore.get_unit('u1')
        assert u1 == unitstore.get_unit('second') * 2

        # Add dimensionless unit
        assert not unitstore.is_defined('u2')
        u2 = unitstore.add_unit('u2', 'dimensionless * 3')
        assert unitstore.is_defined('u2')
        assert str(u2.dimensionality) == 'dimensionless'

        # Make sure 1e6 doesn't get a prefix in it
        unitstore.add_unit('Ms', 'second * 1e6')
        unitstore.add_unit('ks', 'second * 1.e3')

        # Duplicate unit definition
        with pytest.raises(ValueError):
            unitstore.add_unit('u1', 'second * 2')

        # CellML unit redefinition
        with pytest.raises(ValueError):
            unitstore.add_unit('second', 'second * 2')

        # SI prefixes are not recognised
        with pytest.raises(pint.errors.UndefinedUnitError):
            unitstore.add_unit('um', 'micrometer')

        # Pluralising is not switched on
        with pytest.raises(pint.errors.UndefinedUnitError):
            unitstore.add_unit('ms', 'meters')

    def test_add_base_unit(self):
        """Tests UnitStore.add_base_unit()."""

        unitstore = UnitStore()
        assert not unitstore.is_defined('uu')
        uu = unitstore.add_base_unit('uu')
        assert unitstore.is_defined('uu')
        assert unitstore.format(uu, True) in ('1 uu', '1.0 uu')

        # Duplicate unit definition
        with pytest.raises(ValueError):
            unitstore.add_base_unit('uu')

        # CellML unit redefinition
        with pytest.raises(ValueError):
            unitstore.add_base_unit('second')

    def test_convert(self):
        """Tests UnitStore.convert()."""

        store = UnitStore()
        mm = store.add_unit('mm', 'meter / 1000')
        x = 5 * store.get_unit('meter')
        y = store.convert(x, mm)
        assert y == 5000 * mm

    def test_format(self):
        """Tests formatting of units."""

        # Test basic formatting and base unit expansion
        store = UnitStore()
        store.add_unit('x', 'second')
        y = store.add_base_unit('y')
        z = store.add_unit('z', 'y * meter / x ** 2')
        assert store.format(z) == 'z'
        assert store.format(z, True) == '1.0 meter * y / second ** 2'

        # Test if internal name manipulation can mess things up
        a = store.add_base_unit(str(y))
        b = store.add_unit(str(z), '2 * z /' + str(y))
        assert store.format(a) == str(y)
        assert store.format(b) == str(z)
        assert store.format(a, True) in ('1.0 ' + str(y), '1 ' + str(y))
        assert store.format(b, True) == '2.0 meter * y / second ** 2 / ' + str(y)

    def test_get_unit(self):
        """Tests UnitStore.get_unit()."""

        # Get CellML unit
        store = UnitStore()
        assert str(store.get_unit('liter')) == 'liter'
        assert isinstance(store.get_unit('ampere'), store.Unit)

        # Get user unit
        x = store.add_unit('x', 'meter / second')
        assert str(x == 'x')

        # Non-existent unit
        with pytest.raises(KeyError, match='Unknown unit'):
            store.get_unit('towel')

    def test_is_equivalent(self):
        """Tests UnitStore.is_equivalent(a, b)."""

        store = UnitStore()
        a = store.get_unit('meter')
        b = store.add_unit('b', 'meter')
        assert a != b
        assert store.is_equivalent(a, b)

    def test_get_conversion_factor(self):
        """Tests Units.get_conversion_factor() function."""

        store = UnitStore()
        s = store.get_unit('second')
        ms = store.add_unit('ms', 'second / 1000')
        assert store.get_conversion_factor(s, ms) == 0.001
        assert store.get_conversion_factor(ms, s) == 1000

    def test_prefixing(self):
        """Tests edge cases for the internal prefixing."""

        store = UnitStore()
        x = store.add_unit('x', 'meter')

        # Use str(x) as a new unit name: this will have any prefix/postfix the unit store adds so could lead to
        # confusion.
        y = store.add_unit(str(x), 'second')
        z = store.add_unit(str(y), 'ampere')

        # Test use in get_unit()
        assert store.get_unit('x') == x
        assert store.get_unit(str(x)) == y
        assert store.get_unit(str(y)) == z
        assert store.get_unit('x') != store.get_unit(str(x))
        assert store.get_unit('x') != store.get_unit(str(y))
        assert store.get_unit(str(x)) != store.get_unit(str(y))
        assert not store.is_equivalent(store.get_unit('x'), store.get_unit(str(x)))
        assert not store.is_equivalent(store.get_unit('x'), store.get_unit(str(y)))
        assert not store.is_equivalent(store.get_unit(str(x)), store.get_unit(str(y)))

        # Test use in add_unit() expression
        a = store.add_unit('a', str(x))
        b = store.add_unit('b', str(y))
        assert store.is_equivalent(a, y)
        assert store.is_equivalent(b, z)

        # Check that similar names are handled ok
        store = UnitStore()
        a = store.add_unit('my_unit', 'second')
        store.add_unit('also_my_unit', 'volt')
        store.add_unit('my_unit_2', 'ampere')
        assert store.is_equivalent(a, store.get_unit('my_unit'))
        assert not store.is_equivalent(a, store.get_unit('also_my_unit'))
        assert not store.is_equivalent(a, store.get_unit('my_unit_2'))
        assert not store.is_equivalent(store.get_unit('also_my_unit'), store.get_unit('my_unit_2'))

    def test_shared_registry(self):
        """Tests sharing a unit registry."""

        a = UnitStore()
        b = UnitStore(a)

        x = a.add_unit('x', 'meter')
        y = b.add_unit('y', 'meter * 1e-3')
        with pytest.raises(KeyError):
            a.get_unit('y')
        with pytest.raises(KeyError):
            b.get_unit('x')

        assert a.get_conversion_factor(x, y) == 0.001
        assert b.get_conversion_factor(x, y) == 0.001


class TestEvaluateUnits(object):
    """Tests UnitStore.evaluate_units()."""

    def test_numbers(self):
        """Tests on numbers."""

        store = UnitStore()
        _1 = NumberDummy(1, store.get_unit('kelvin'))
        _2 = NumberDummy(2, store.get_unit('dimensionless'))

        assert store.evaluate_units(_1) == store.get_unit('kelvin')
        assert store.evaluate_units(_2) == store.get_unit('dimensionless')
        assert store.evaluate_units(sp.sympify('1.0')) == store.get_unit('dimensionless')

    # TODO: Split into multiple methods and remove redundant tests
    def test_others(self):
        """ Tests the units.UnitCalculator class. """

        store = UnitStore()
        a = VariableDummy('a', store.get_unit('meter'))
        b = VariableDummy('b', store.get_unit('second'))
        c = VariableDummy('c', store.get_unit('gram'))
        d = VariableDummy('d', store.get_unit('meter'))
        x = VariableDummy('x', store.get_unit('kilogram'))
        y = VariableDummy('y', store.get_unit('volt'))
        n = VariableDummy('n', store.get_unit('dimensionless'))
        av = VariableDummy('av', store.get_unit('meter'), initial_value=2)
        _1 = NumberDummy(1, store.get_unit('kelvin'))
        _2 = NumberDummy(2, store.get_unit('dimensionless'))
        _25 = NumberDummy(2.5, store.get_unit('dimensionless'))

        # functions were all args should have same unit
        # and this is the unit that will be returned
        # pass cases
        assert store.evaluate_units(-a) == store.get_unit('meter')
        assert store.evaluate_units(a + a + a + a) == store.get_unit('meter')
        assert store.evaluate_units(a + 2 * a + 3 * a + 4 * a) == store.get_unit('meter')
        assert store.evaluate_units(2 * a - a) == store.get_unit('meter')
        assert store.evaluate_units(a + d) == store.get_unit('meter')
        assert store.evaluate_units(a - d) == store.get_unit('meter')
        assert store.evaluate_units(sp.Abs(-2 * y)) == store.get_unit('volt')
        assert store.evaluate_units(sp.floor(x)) == store.get_unit('kilogram')
        assert store.evaluate_units(sp.floor(_1)) == store.get_unit('kelvin')
        assert store.evaluate_units(sp.floor(12.5) * a) == store.get_unit('meter')
        assert store.evaluate_units(sp.ceiling(x)) == store.get_unit('kilogram')
        assert store.evaluate_units(sp.ceiling(12.6) * b) == store.get_unit('second')
        assert store.evaluate_units(sp.ceiling(_1)) == store.get_unit('kelvin')

        # functions were all args should have same unit
        # and this is the unit that will be returned
        # fail cases
        with pytest.raises(InputArgumentsInvalidUnitsError):
            store.evaluate_units(a + b)
        try:
            store.evaluate_units(a + b)
        except InputArgumentsInvalidUnitsError as err:
            assert err.message == 'The arguments to this expression should all have the same units.'
            assert err.expression == '_a + _b'
        with pytest.raises(InputArgumentsInvalidUnitsError):
            store.evaluate_units(a - b)
        try:
            store.evaluate_units(a - b)
        except InputArgumentsInvalidUnitsError as err:
            assert err.message == 'The arguments to this expression should all have the same units.'
            assert err.expression == '_a - _b'
        try:
            sp.floor(x, y)
        except TypeError as err:
            assert err.args[0] == 'floor takes exactly 1 argument (2 given)'

        # special case - piecewise
        m = store.evaluate_units(sp.Piecewise((a, x < 1), (a + a, x > 1), (3 * a, True)))
        assert m == store.get_unit('meter')
        # fail special case -piecewise
        with pytest.raises(InputArgumentsInvalidUnitsError):
            store.evaluate_units(sp.Piecewise((a, x < 1), (b, x > 1), (c, True)))
        try:
            store.evaluate_units(sp.Piecewise((a, x < 1), (b, x > 1), (c, True)))
        except InputArgumentsInvalidUnitsError as err:
            assert err.message == 'The arguments to this expression should all have the same units.'
            assert err.expression == 'Piecewise((_a, _x < 1), (_b, _x > 1), (_c, True))'

        # cases with any units allowed as arguments
        assert store.evaluate_units((a * a) / b) == (
            store.get_unit('meter') ** 2 / store.get_unit('second'))

        # root and power
        assert store.evaluate_units(a**_2) == store.get_unit('meter')**2
        assert store.evaluate_units(sp.sqrt(c ** 2)) == store.get_unit('gram')
        assert store.evaluate_units(sp.sqrt(a * d)) == store.get_unit('meter')
        # root and power fails
        expr = a ** _1
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            store.evaluate_units(expr)
        try:
            store.evaluate_units(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The second argument to this expression should be dimensionless.'
            assert err.expression == '_a**_1'
        expr = a + a ** n
        with pytest.raises(InputArgumentMustBeNumberError):
            store.evaluate_units(expr)
        try:
            store.evaluate_units(expr)
        except InputArgumentMustBeNumberError as err:
            assert err.message == 'The second argument to this expression should be a number.'
            assert err.expression == '_a**_n'

        # bizarre cases
        expr = av ** _25
        result = store.evaluate_units(expr)
        assert result == store.get_unit('meter')**2.5  # and result.magnitude == 5.6568
        expr = sp.root(av, _25)
        result = store.evaluate_units(expr)
        assert result == store.get_unit('meter')**0.4  # and result.magnitude == 1.319507

        # relational operators throw an exception
        expr = a > _1
        with pytest.raises(BooleanUnitsError):
            store.evaluate_units(expr)

        # special case - derivative
        dadb = sp.diff(a * b, b)
        assert store.evaluate_units(dadb) == store.get_unit('meter')

        dadb = sp.Derivative(a * b, b)
        assert store.evaluate_units(dadb) == store.get_unit('meter')

        # log and exponent
        assert store.evaluate_units(sp.log(_2)) == store.get_unit('dimensionless')
        assert store.evaluate_units(sp.log(n)) == store.get_unit('dimensionless')
        assert store.evaluate_units(sp.ln(_2)) == store.get_unit('dimensionless')
        assert store.evaluate_units(sp.log(n, _2)) == store.get_unit('dimensionless')
        assert store.evaluate_units(sp.log(_2, 3)) == store.get_unit('dimensionless')
        assert store.evaluate_units(sp.exp(_2)) == store.get_unit('dimensionless')
        assert store.evaluate_units(sp.exp(n)) == store.get_unit('dimensionless')
        # log and exponent - fails
        expr = sp.log(a)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            store.evaluate_units(expr)
        try:
            store.evaluate_units(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The arguments to this expression should be dimensionless.'
            assert err.expression == 'log(_a)'
        expr = sp.log(_2, av)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            store.evaluate_units(expr)
        try:
            store.evaluate_units(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The arguments to this expression should be dimensionless.'
            assert err.expression == 'log(_av)'
        expr = sp.exp(av)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            store.evaluate_units(expr)
        try:
            store.evaluate_units(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The arguments to this expression should be dimensionless.'
            assert err.expression == 'exp(_av)'

        # trig functions
        assert store.evaluate_units(sp.sin(n)) == store.get_unit('dimensionless')
        # trig functions - fails
        expr = sp.cos(av)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            store.evaluate_units(expr)
        try:
            store.evaluate_units(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The arguments to this expression should be dimensionless.'
            assert err.expression == 'cos(_av)'

        # constants
        assert store.evaluate_units(sp.pi) == store.get_unit('dimensionless')
        assert store.evaluate_units(sp.E) == store.get_unit('dimensionless')
        assert store.evaluate_units(sp.oo) == store.get_unit('dimensionless')
        assert store.evaluate_units(sp.nan) == store.get_unit('dimensionless')

        expr = sp.true
        with pytest.raises(BooleanUnitsError):
            store.evaluate_units(expr)
        try:
            store.evaluate_units(expr)
        except BooleanUnitsError as err:
            assert err.message == 'This expression involves boolean values which do not conform to ' \
                                  'unit dimensionality rules.'
            assert err.expression == 'True'

        # factorial needs its own catch
        assert store.evaluate_units(sp.factorial(n)) == store.get_unit('dimensionless')
        expr = sp.factorial(av)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            store.evaluate_units(expr)
        try:
            store.evaluate_units(expr)
        except InputArgumentsMustBeDimensionlessError as err:
            assert err.message == 'The arguments to this expression should be dimensionless.'
            assert err.expression == 'factorial(_av)'

        # logic functions throw an exception
        expr = a & b
        with pytest.raises(BooleanUnitsError):
            store.evaluate_units(expr)
        try:
            store.evaluate_units(expr)
        except BooleanUnitsError as err:
            assert err.message == 'This expression involves boolean values which do not conform to ' \
                                  'unit dimensionality rules.'
            assert err.expression == '_a & _b'

        # check that not gets caught as it is listed separately
        # does indeed get caught by is_Boolean
        expr = ~a
        with pytest.raises(BooleanUnitsError):
            store.evaluate_units(expr)

        # this is a generic test of a function with one dimensionless argument
        # this uses a function that is not in cellml - just to check the logic
        assert store.evaluate_units(sp.sign(n)) == store.get_unit('dimensionless')

        expr = sp.Matrix([[1, 0], [0, 1]])
        with pytest.raises(UnexpectedMathUnitsError):
            store.evaluate_units(expr)

        N = sp.Matrix([[1, 0], [0, 1]])
        M = sp.Matrix([[1, 0], [0, 1]])
        expr = M + N
        with pytest.raises(UnexpectedMathUnitsError):
            store.evaluate_units(expr)

        expr = sp.cos(x).series(x, 0, 10)
        with pytest.raises(UnexpectedMathUnitsError):
            store.evaluate_units(expr)
