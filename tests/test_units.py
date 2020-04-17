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
    UnitConversionError,
    UnitStore,
)


@pytest.fixture(scope='class')
def store():
    return UnitStore()


@pytest.fixture(scope='class')
def cm(store):
    return store.add_unit('cm', '0.01 * metre')


@pytest.fixture(scope='class')
def ms(store):
    return store.add_unit('ms', '0.001 * second')


@pytest.fixture(scope='class')
def microsecond(store):
    return store.add_unit('microsecond', '1e-6 * second')


class TestUnits:

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


class TestEvaluateUnits:
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


class TestConvertingExpressions:
    """Test the UnitStore.convert_expression_recursively and set_lhs_units_from_rhs methods."""
    def test_variable_no_conversion(self, store):
        x = VariableDummy('x', store.get_unit('metre'))

        new_x = store.convert_expression_recursively(x, None)
        assert x is new_x

        new_x = store.convert_expression_recursively(x, x.units)
        assert x is new_x

    def test_variable_conversion(self, store, cm):
        x = VariableDummy('x', store.get_unit('metre'))
        new_x = store.convert_expression_recursively(x, cm)
        assert str(new_x) == '_100.0*_x'
        assert new_x.args[1] is x
        assert isinstance(new_x.args[0], NumberDummy)
        assert new_x.args[0].units == cm / store.get_unit('metre')

    def test_number_conversion(self, store, cm):
        _5 = NumberDummy(5, cm)
        new_5 = store.convert_expression_recursively(_5, store.get_unit('metre'))
        assert str(new_5) == '_0.01*_5'
        assert new_5.args[1] is _5
        assert isinstance(new_5.args[0], NumberDummy)
        assert new_5.args[0].units == store.get_unit('metre') / cm

    def test_plain_numbers(self, store):
        dimensionless = store.get_unit('dimensionless')
        assert store.convert_expression_recursively(sp.E, dimensionless) is sp.E
        assert store.convert_expression_recursively(sp.pi, None) is sp.pi
        assert store.convert_expression_recursively(sp.oo, dimensionless) is sp.oo
        assert store.convert_expression_recursively(sp.nan, None) is sp.nan
        assert store.convert_expression_recursively(sp.true, dimensionless) is sp.true
        assert store.convert_expression_recursively(sp.false, None) is sp.false

        expr = sp.Integer(2)
        assert store.convert_expression_recursively(expr, dimensionless) is expr

        expr = sp.Rational(2, 3)
        assert store.convert_expression_recursively(expr, None) is expr

    def test_derivative_no_conversion(self, store):
        x = VariableDummy('x', store.get_unit('metre'))
        t = VariableDummy('t', store.get_unit('second'))
        expr = sp.Derivative(x, t)

        new_expr = store.convert_expression_recursively(expr, None)
        assert expr is new_expr

        new_expr = store.convert_expression_recursively(expr, x.units / t.units)
        assert expr is new_expr

    def test_derivative_conversion(self, store, cm, ms):
        x = VariableDummy('x', store.get_unit('metre'))
        t = VariableDummy('t', store.get_unit('second'))
        expr = sp.Derivative(x, t)

        new_expr = store.convert_expression_recursively(expr, cm / t.units)
        assert str(new_expr) == '_100.0*Derivative(_x, _t)'
        assert new_expr.args[0].args[0] is x
        assert new_expr.args[0].args[1][0] is t

        new_expr = store.convert_expression_recursively(expr, x.units / ms)
        assert str(new_expr) == '_0.001*Derivative(_x, _t)'

    def test_mul_and_pow(self, store, cm, ms):
        x = VariableDummy('x', cm)
        y = VariableDummy('y', ms)
        expr = x / y  # Becomes x * (1/y)

        # No conversion
        new_expr = store.convert_expression_recursively(expr, None)
        assert expr is new_expr

        # With conversion
        new_expr = store.convert_expression_recursively(expr, store.get_unit('metre') / store.get_unit('second'))
        assert str(new_expr) == '_10.0*_x/_y'
        assert new_expr.args[2] is x
        assert new_expr.args[0].args[0] is y

        # With conversion only for exponent
        _4 = NumberDummy('4', store.get_unit('second'))
        _2 = NumberDummy('2000', ms)
        expr = x ** (_4 / _2)
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == '_x**(_1000.0*_4/_2000)'

        # With a base that needs internal conversion
        expr = (y + _4) ** 2
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == '(_0.001*_y + _4)**2'

    def test_square_root(self, store, ms):
        x = VariableDummy('x', store.get_unit('second') ** 2)
        expr = x ** (1 / 2)
        new_expr = store.convert_expression_recursively(expr, ms)
        assert str(new_expr) == '_1000.0*_x**0.5'
        assert new_expr.args[0].args[0] is x

    def test_add_and_subtract(self, store, ms, microsecond):
        x = VariableDummy('x', store.get_unit('second'))
        y = VariableDummy('y', ms)
        z = VariableDummy('z', microsecond)

        # If no conversion is needed we get the original expression
        expr = y + NumberDummy('2', ms)
        assert store.convert_expression_recursively(expr, ms) is expr

        # If we don't specify units, the first argument (y in canonical form) is chosen
        expr = z + x - y
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == '_0.001*_z + _1000.0*_x - _y'

        new_expr = store.convert_expression_recursively(expr, microsecond)
        assert str(new_expr) == '-_1000.0*_y + _1000000.0*_x + _z'

    def test_abs_ceil_floor(self, store, ms):
        x = VariableDummy('x', store.get_unit('second'))

        expr = sp.Abs(x)
        new_expr = store.convert_expression_recursively(expr, None)
        assert new_expr is expr

        expr = sp.Abs(x)
        new_expr = store.convert_expression_recursively(expr, ms)
        assert isinstance(new_expr, sp.Abs)
        assert str(new_expr) == 'Abs(_1000.0*_x)'
        assert new_expr.args[0].args[1] is x

        expr = sp.floor(x)
        new_expr = store.convert_expression_recursively(expr, store.get_unit('second'))
        assert new_expr is expr

        expr = sp.floor(x)
        new_expr = store.convert_expression_recursively(expr, ms)
        assert isinstance(new_expr, sp.floor)
        assert str(new_expr) == 'floor(_1000.0*_x)'
        assert new_expr.args[0].args[1] is x

        expr = sp.ceiling(x)
        new_expr = store.convert_expression_recursively(expr, ms)
        assert isinstance(new_expr, sp.ceiling)
        assert str(new_expr) == 'ceiling(_1000.0*_x)'
        assert new_expr.args[0].args[1] is x

    def test_exp_log_trig(self, store, ms):
        dimensionless = store.get_unit('dimensionless')
        x = VariableDummy('x', dimensionless)
        y = VariableDummy('y', store.get_unit('second'))
        z = VariableDummy('z', ms)

        expr = sp.exp(x)
        assert store.convert_expression_recursively(expr, dimensionless) is expr

        expr = sp.log(y / z)
        new_expr = store.convert_expression_recursively(expr, None)
        assert isinstance(new_expr, sp.log)
        assert str(new_expr) == 'log(_1000.0*_y/_z)'

        expr = sp.sin(z / y)
        new_expr = store.convert_expression_recursively(expr, dimensionless)
        assert isinstance(new_expr, sp.sin)
        assert str(new_expr) == 'sin(_0.001*_z/_y)'

    def test_relations(self, store, ms):
        x = VariableDummy('x', ms)
        y = VariableDummy('y', store.get_unit('second'))
        z = VariableDummy('z', ms)

        expr = y < z
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == '_y < _0.001*_z'

        expr = z > y
        new_expr = store.convert_expression_recursively(expr, store.get_unit('dimensionless'))
        assert str(new_expr) == '_z > _1000.0*_y'

        # Case with no conversion
        expr = x < z
        new_expr = store.convert_expression_recursively(expr, None)
        assert new_expr is expr

    def test_piecewise(self, store, ms, microsecond, cm):
        # This also checks more complex nested expressions
        dimensionless = store.get_unit('dimensionless')
        _2 = NumberDummy(2, dimensionless)
        x = VariableDummy('x', dimensionless)
        y = VariableDummy('y', store.get_unit('second'))
        z = VariableDummy('z', ms)

        expr = sp.Piecewise(
            (y, sp.And(y < z, x > _2)),
            (z, x < _2),
            (z * _2, True),
        )

        # Units of result will be chosen from first case, i.e. second
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == (
            'Piecewise((_y, (_2 < _x) & (_y < _0.001*_z)), (_0.001*_z, _2 > _x), (_0.001*_2*_z, True))')

        # Units of result are specified
        new_expr = store.convert_expression_recursively(expr, ms)
        assert str(new_expr) == (
            'Piecewise((_1000.0*_y, (_2 < _x) & (_y < _0.001*_z)), (_z, _2 > _x), (_2*_z, True))')

        # A simpler case with no conversion
        _1 = NumberDummy('1', ms)
        expr = sp.Piecewise((z, x < _2), (_1, True))
        assert store.convert_expression_recursively(expr, ms) is expr

    def test_assignment(self, store, ms):
        x = VariableDummy('x', store.get_unit('second'))
        _10 = NumberDummy(10, ms)
        expr = sp.Eq(x, _10)
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == 'Eq(_x, _0.001*_10)'

    def test_set_lhs_units_from_rhs(self, store, ms):
        x = VariableDummy('x', None)
        _10 = NumberDummy(10, ms)
        expr = sp.Eq(x, _10)
        new_expr = store.set_lhs_units_from_rhs(expr)
        assert str(new_expr) == 'Eq(_x, _10)'
        assert new_expr.args[0] is x
        assert new_expr.args[1] is _10
        assert x.units is ms

    def test_set_lhs_units_from_converted_rhs(self, store, ms):
        x = VariableDummy('x', None)
        y = VariableDummy('y', store.get_unit('second'))
        _10 = NumberDummy(10, ms)
        expr = sp.Eq(x, _10 + y)
        new_expr = store.set_lhs_units_from_rhs(expr)
        assert str(new_expr) == 'Eq(_x, _10 + _1000.0*_y)'
        assert new_expr.args[0] is x
        assert x.units is ms

    # Methods below check error cases
    def test_symbol_wrong_dimensions(self, store):
        x = VariableDummy('x', store.get_unit('metre'))
        with pytest.raises(UnitConversionError, match='from meter to second'):
            store.convert_expression_recursively(x, store.get_unit('second'))

    def test_derivative_wrong_dimensions(self, store):
        x = VariableDummy('x', store.get_unit('metre'))
        t = VariableDummy('t', store.get_unit('second'))
        expr = sp.Derivative(x, t)
        with pytest.raises(UnitConversionError, match='Context: trying to convert'):
            store.convert_expression_recursively(expr, store.get_unit('metre'))

    def test_complex_derivative(self, store):
        x = VariableDummy('x', store.get_unit('metre'))
        y = VariableDummy('y', store.get_unit('metre'))
        t = VariableDummy('t', store.get_unit('second'))
        t2 = VariableDummy('t2', store.get_unit('second'))
        expr = sp.Derivative(x + y, t, t2)
        with pytest.raises(UnexpectedMathUnitsError,
                           match='only support first order derivatives of single variables'):
            store.convert_expression_recursively(expr, store.get_unit('metre'))

    def test_mul_wrong_dimensions(self, store):
        x = VariableDummy('x', store.get_unit('metre'))
        _1 = NumberDummy(1, store.get_unit('second'))
        expr = x * _1
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, store.get_unit('metre'))

    def test_error_matrix(self, store):
        expr = sp.Matrix([[1, 0], [0, 1]])
        with pytest.raises(UnexpectedMathUnitsError):
            store.convert_expression_recursively(expr, None)

    def test_error_exponent_not_dimensionless(self, store):
        x = VariableDummy('x', store.get_unit('metre'))
        _1 = NumberDummy(1, store.get_unit('second'))
        expr = x ** _1
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, None)

    def test_error_exponent_not_number(self, store):
        x = VariableDummy('x', store.get_unit('dimensionless'))
        _1 = NumberDummy(1, store.get_unit('second'))
        expr = _1 ** x
        with pytest.raises(InputArgumentMustBeNumberError):
            store.convert_expression_recursively(expr, None)

    def test_pow_wrong_dimensions(self, store):
        x = VariableDummy('x', store.get_unit('second'))
        _2 = NumberDummy(2, store.get_unit('dimensionless'))
        expr = x ** _2
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, store.get_unit('second'))

    def test_add_wrong_dimensions(self, store):
        x = VariableDummy('x', store.get_unit('second'))
        _2 = NumberDummy(2, store.get_unit('dimensionless'))
        expr = x + _2
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, None)

    def test_relational_must_be_dimensionless(self, store):
        x = VariableDummy('x', store.get_unit('second'))
        y = VariableDummy('y', store.get_unit('second'))
        expr = x < y
        with pytest.raises(BooleanUnitsError):
            store.convert_expression_recursively(expr, store.get_unit('second'))

    def test_relational_dimension_mismatch(self, store):
        x = VariableDummy('x', store.get_unit('second'))
        y = VariableDummy('y', store.get_unit('metre'))
        expr = x <= y
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, store.get_unit('dimensionless'))

    def test_piecewise_condition_not_dimensionless(self, store):
        a = VariableDummy('a', store.get_unit('metre'))
        x = VariableDummy('x', store.get_unit('second'))
        expr = sp.Piecewise((a, x), (-a, True))
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, store.get_unit('metre'))

    def test_piecewise_cannot_convert_result(self, store):
        a = VariableDummy('a', store.get_unit('metre'))
        x = VariableDummy('x', store.get_unit('dimensionless'))
        expr = sp.Piecewise((a, x), (-a, True))
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, store.get_unit('second'))

    def test_exp_must_be_dimensionless(self, store):
        x = VariableDummy('x', store.get_unit('dimensionless'))
        expr = sp.exp(x)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            store.convert_expression_recursively(expr, store.get_unit('second'))

    def test_number_must_be_dimensionless(self, store):
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            store.convert_expression_recursively(sp.E, store.get_unit('second'))

    def test_error_unsupported_math(self, store):
        x = VariableDummy('x', store.get_unit('second'))
        expr = sp.Integral(x**2, x)
        with pytest.raises(UnexpectedMathUnitsError):
            store.convert_expression_recursively(expr, None)

    def test_set_lhs_units_from_inconsistent_rhs(self, store):
        x = VariableDummy('x', None)
        _10 = NumberDummy(10, store.get_unit('metre'))
        _20 = NumberDummy(20, store.get_unit('second'))
        expr = sp.Eq(x, _10 + _20)
        with pytest.raises(UnitConversionError, match='Context: trying to set LHS units'):
            new_expr = store.set_lhs_units_from_rhs(expr)
