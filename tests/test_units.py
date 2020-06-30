import pint
import pytest
import sympy as sp

from cellmlmanip.model import Quantity, Variable
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


@pytest.fixture(scope='class')
def a(store):
    return Variable('a', store.get_unit('meter'))


@pytest.fixture(scope='class')
def b(store):
    return Variable('b', store.get_unit('second'))


@pytest.fixture(scope='class')
def c(store):
    return Variable('c', store.get_unit('gram'))


@pytest.fixture(scope='class')
def d(store):
    return Variable('d', store.get_unit('meter'))


@pytest.fixture(scope='class')
def x(store):
    return Variable('x', store.get_unit('kilogram'))


@pytest.fixture(scope='class')
def y(store):
    return Variable('y', store.get_unit('volt'))


@pytest.fixture(scope='class')
def n(store):
    return Variable('n', store.get_unit('dimensionless'))


@pytest.fixture(scope='class')
def av(store):
    return Variable('av', store.get_unit('meter'), initial_value=2)


@pytest.fixture(scope='class')
def _1(store):
    return Quantity(1, store.get_unit('kelvin'))


@pytest.fixture(scope='class')
def _2(store):
    return Quantity(2, store.get_unit('dimensionless'))


@pytest.fixture(scope='class')
def _25(store):
    return Quantity(2.5, store.get_unit('dimensionless'))


@pytest.fixture(scope='class')
def mtr(store):
    return Variable('mtr', store.get_unit('metre'))


@pytest.fixture(scope='class')
def sec(store):
    return Variable('sec', store.get_unit('second'))


@pytest.fixture(scope='class')
def dim(store):
    return Variable('dim', store.get_unit('dimensionless'))


class TestUnits:

    def test_add_unit(self, store):
        """Tests UnitStore.add_unit()."""

        # Add ordinary unit
        assert not store.is_defined('u1')
        u1 = store.add_unit('u1', 'second * 2')
        assert store.is_defined('u1')
        assert u1 == store.get_unit('u1')
        assert u1 == store.get_unit('second') * 2

        # Add dimensionless unit
        assert not store.is_defined('u2')
        u2 = store.add_unit('u2', 'dimensionless * 3')
        assert store.is_defined('u2')
        assert str(u2.dimensionality) == 'dimensionless'

        # Make sure 1e6 doesn't get a prefix in it
        store.add_unit('Ms', 'second * 1e6')
        store.add_unit('ks', 'second * 1.e3')

        # Duplicate unit definition
        with pytest.raises(ValueError):
            store.add_unit('u1', 'second * 2')

        # CellML unit redefinition
        with pytest.raises(ValueError):
            store.add_unit('second', 'second * 2')

        # SI prefixes are not recognised
        with pytest.raises(pint.errors.UndefinedUnitError):
            store.add_unit('um', 'micrometer')

        # Pluralising is not switched on
        with pytest.raises(pint.errors.UndefinedUnitError):
            store.add_unit('ms', 'meters')

    def test_add_base_unit(self, store):
        """Tests UnitStore.add_base_unit()."""

        assert not store.is_defined('uu')
        uu = store.add_base_unit('uu')
        assert store.is_defined('uu')
        assert store.format(uu, True) in ('1 uu', '1.0 uu')

        # Duplicate unit definition
        with pytest.raises(ValueError):
            store.add_base_unit('uu')

        # CellML unit redefinition
        with pytest.raises(ValueError):
            store.add_base_unit('second')

    def test_convert(self, store):
        """Tests UnitStore.convert()."""

        mm = store.add_unit('mm', 'meter / 1000')
        qx = 5 * store.get_unit('meter')
        qy = store.convert(qx, mm)
        assert qy == 5000 * mm

    def test_format(self, store):
        """Tests formatting of units."""

        # Test basic formatting and base unit expansion
        store.add_unit('se', 'second')
        y = store.add_base_unit('y')
        z = store.add_unit('z', 'y * meter / se ** 2')
        assert store.format(z) == 'z'
        assert store.format(z, True) == '1.0 meter * y / second ** 2'

        # Test if internal name manipulation can mess things up
        a = store.add_base_unit(str(y))
        b = store.add_unit(str(z), '2 * z /' + str(y))
        assert store.format(a) == str(y)
        assert store.format(b) == str(z)
        assert store.format(a, True) in ('1.0 ' + str(y), '1 ' + str(y))
        assert store.format(b, True) == '2.0 meter * y / second ** 2 / ' + str(y)

    def test_get_unit(self, store):
        """Tests UnitStore.get_unit()."""

        # Get CellML unit
        assert str(store.get_unit('liter')) == 'liter'
        assert isinstance(store.get_unit('ampere'), store.Unit)

        # Get user unit
        m_per_s = store.add_unit('m_per_s', 'meter / second')
        get_m_per_s = store.get_unit('m_per_s')
        assert str(m_per_s == get_m_per_s) and str(get_m_per_s).endswith('m_per_s')

        # Non-existent unit
        with pytest.raises(KeyError, match='Unknown unit'):
            store.get_unit('towel')

    def test_is_equivalent(self, store):
        """Tests UnitStore.is_equivalent(a, b)."""

        a = store.get_unit('meter')
        b = store.add_unit('b', 'meter')
        assert a != b
        assert store.is_equivalent(a, b)

    def test_get_conversion_factor(self, store):
        """Tests Units.get_conversion_factor() function."""

        s = store.get_unit('second')
        ms = store.add_unit('ms', 'second / 1000')
        assert store.get_conversion_factor(ms, s) == 0.001
        assert store.get_conversion_factor(s, ms) == 1000

    def test_prefixing(self, store):
        """Tests edge cases for the internal prefixing."""

        mtr2 = store.add_unit('mtr2', 'meter')

        # Use str(x) as a new unit name: this will have any prefix/postfix the unit store adds so could lead to
        # confusion.
        y = store.add_unit(str(mtr2), 'second')
        z = store.add_unit(str(y), 'ampere')

        # Test use in get_unit()
        assert store.get_unit('mtr2') == mtr2
        assert store.get_unit(str(mtr2)) == y
        assert store.get_unit(str(y)) == z
        assert store.get_unit('mtr2') != store.get_unit(str(mtr2))
        assert store.get_unit('mtr2') != store.get_unit(str(y))
        assert store.get_unit(str(mtr2)) != store.get_unit(str(y))
        assert not store.is_equivalent(store.get_unit('mtr2'), store.get_unit(str(mtr2)))
        assert not store.is_equivalent(store.get_unit('mtr2'), store.get_unit(str(y)))
        assert not store.is_equivalent(store.get_unit(str(mtr2)), store.get_unit(str(y)))

        # Test use in add_unit() expression
        a = store.add_unit('a', str(mtr2))
        b2 = store.add_unit('b2', str(y))
        assert store.is_equivalent(a, y)
        assert store.is_equivalent(b2, z)

        # Check that similar names are handled ok
        store = UnitStore()
        a = store.add_unit('my_unit', 'second')
        store.add_unit('also_my_unit', 'volt')
        store.add_unit('my_unit_2', 'ampere')
        assert store.is_equivalent(a, store.get_unit('my_unit'))
        assert not store.is_equivalent(a, store.get_unit('also_my_unit'))
        assert not store.is_equivalent(a, store.get_unit('my_unit_2'))
        assert not store.is_equivalent(store.get_unit('also_my_unit'), store.get_unit('my_unit_2'))

    def test_shared_registry(self, store, mtr):
        """Tests sharing a unit registry."""

        a = store
        b = UnitStore(a)

        mtr = mtr.units
        yy = b.add_unit('yy', 'meter * 1e-3')
        with pytest.raises(KeyError):
            a.get_unit('yy')
        with pytest.raises(KeyError):
            b.get_unit('mtr')

        assert a.get_conversion_factor(mtr, yy) == 1000
        assert b.get_conversion_factor(mtr, yy) == 1000


class TestEvaluateUnits:
    """Tests UnitStore.evaluate_units()."""

    def test_numbers(self, store, _1, _2):
        """Tests on numbers."""

        assert store.evaluate_units(_1) == store.get_unit('kelvin')
        assert store.evaluate_units(_2) == store.get_unit('dimensionless')
        assert store.evaluate_units(sp.sympify('1.0')) == store.get_unit('dimensionless')

    def test_func_args_same_units_(self, store, a, b, d, x, y, _1):
        """ Tests functions were all args should have same unit and this is the unit returned. """

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

    def test_func_args_same_units_fail_cases(self, store, a, b, x, y):
        """ Tests functions were all args should have same unit and this is the unit returned, fail cases. """

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

    def test_special_piecewise(self, store, a, b, c, x):
        """ Tests special case - piecewise """

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

    def test_any_units_as_args(self, store, a, b):
        """ Tests cases with any units allowed as arguments. """

        assert store.evaluate_units((a * a) / b) == (
            store.get_unit('meter') ** 2 / store.get_unit('second'))

    def test_root_and_power(self, store, a, c, d, n, av, _1, _2, _25):
        """ Tests root and power. """

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

        # Fractional powers via quantities
        expr = av ** _25
        result = store.evaluate_units(expr)
        assert result == store.get_unit('meter')**2.5
        expr = sp.root(av, _25)
        result = store.evaluate_units(expr)
        assert result == store.get_unit('meter')**0.4

    def test_relational_operators(self, caplog, store, a, d, _1):
        """ Tests relational operators throw an exception, unequal units of quantities triggers a log message. """

        expr = a > _1
        with pytest.raises(BooleanUnitsError):
            store.evaluate_units(expr)
        assert 'Relational args do not have the same unit: _a > _1' in caplog.text

        # relational operators throw an exception
        # units of quantities are equal
        expr = a > d
        with pytest.raises(BooleanUnitsError):
            store.evaluate_units(expr)
        assert 'Relational args do not have the same unit: _a > _s' not in caplog.text

    def test_derivative(self, store, a, b, c, d, x, y, n, av, _1, _2, _25):
        """ Tests special case - derivative. """

        dadb = sp.diff(a * b, b)
        assert store.evaluate_units(dadb) == store.get_unit('meter')

        dadb = sp.Derivative(a * b, b)
        assert store.evaluate_units(dadb) == store.get_unit('meter')

    def test_log_and_exponent(self, store, a, n, av, _2):
        """ Tests log and exponent. """

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

    def test_trig_functions(self, store, n, av):
        """ Tests trig functions. """

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

    def test_constants(self, store):
        """ Tests constants. """

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

    def test_factorial(self, store, n, av):
        """ Tests factorial. """

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

    def test_logic_functions(self, store, a, b):
        """ Tests logic functions. """

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

    def test_function_diemensionless_arg(self, store, x, n):
        """ This is a generic test of a function with one dimensionless argument
            this uses a function that is not in cellml - just to check the logic. """

        assert store.evaluate_units(sp.sign(n)) == store.get_unit('dimensionless')

        # Check we get an error if evaluating units for function with multiple dimensionless arguments
        with pytest.raises(UnexpectedMathUnitsError):
            assert store.evaluate_units(sp.Function('sgn')(n, n)) == store.get_unit('dimensionless')

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
    """
    Test the UnitStore methods that perform within-equation conversion.
    """
    def test_variable_no_conversion(self, store, mtr):

        new_mtr = store.convert_expression_recursively(mtr, None)
        assert mtr is new_mtr

        new_mtr = store.convert_expression_recursively(mtr, mtr.units)
        assert mtr is new_mtr

    def test_variable_conversion(self, store, mtr, cm):
        new_mtr = store.convert_expression_recursively(mtr, cm)
        assert str(new_mtr) == '_100*_mtr'
        assert new_mtr.args[1] is mtr
        assert isinstance(new_mtr.args[0], Quantity)
        assert new_mtr.args[0].units == cm / store.get_unit('metre')

    def test_number_conversion(self, store, cm):
        _5 = Quantity(5, cm)
        new_5 = store.convert_expression_recursively(_5, store.get_unit('metre'))
        assert str(new_5) == '_0.01*_5'
        assert new_5.args[1] is _5
        assert isinstance(new_5.args[0], Quantity)
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

    def test_derivative_no_conversion(self, store, mtr):
        t = Variable('t', store.get_unit('second'))
        expr = sp.Derivative(mtr, t)

        new_expr = store.convert_expression_recursively(expr, None)
        assert expr is new_expr

        new_expr = store.convert_expression_recursively(expr, mtr.units / t.units)
        assert expr is new_expr

    def test_derivative_conversion(self, store, mtr, cm, ms):
        t = Variable('t', store.get_unit('second'))
        expr = sp.Derivative(mtr, t)

        new_expr = store.convert_expression_recursively(expr, cm / t.units)
        assert str(new_expr) == '_100*Derivative(_mtr, _t)'
        assert new_expr.args[0].args[0] is mtr
        assert new_expr.args[0].args[1][0] is t

        new_expr = store.convert_expression_recursively(expr, mtr.units / ms)
        assert str(new_expr) == '_0.001*Derivative(_mtr, _t)'

    def test_mul_and_pow(self, store, cm, ms):
        cmtr = Variable('cmtr', cm)
        y = Variable('y', ms)
        expr = cmtr / y  # Becomes cmtr * (1/y)

        # No conversion
        new_expr = store.convert_expression_recursively(expr, None)
        assert expr is new_expr

        # With conversion
        new_expr = store.convert_expression_recursively(expr, store.get_unit('metre') / store.get_unit('second'))
        assert str(new_expr) == '_10*_cmtr/_y'
        assert new_expr.args[2] is cmtr
        assert new_expr.args[0].args[0] is y

        # With conversion only for exponent
        _4 = Quantity('4', store.get_unit('second'))
        _2 = Quantity('2000', ms)
        expr = cmtr ** (_4 / _2)
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == '_cmtr**(_1000*_4/_2000)'

        # With a base that needs internal conversion
        expr = (y + _4) ** 2
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == '(_0.001*_y + _4)**2'

    def test_mul_with_converted_arg(self, store, sec, ms):
        y = Variable('y', ms)
        z = Variable('z', ms)
        expr = (sec + y) * z
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == '_z*(_0.001*_y + _sec)'

    def test_square_root(self, store, ms):
        secsqr = Variable('secsqr', store.get_unit('second') ** 2)
        expr = secsqr ** (1 / 2)
        new_expr = store.convert_expression_recursively(expr, ms)
        assert str(new_expr) == '_1000*_secsqr**0.5'
        assert new_expr.args[0].args[0] is secsqr

    def test_add_and_subtract(self, store, sec, ms, microsecond):
        y = Variable('y', ms)
        z = Variable('z', microsecond)

        # If no conversion is needed we get the original expression
        expr = y + Quantity('2', ms)
        assert store.convert_expression_recursively(expr, ms) is expr

        # If we don't specify units, the first argument (y in canonical form) is chosen
        expr = z + sec - y
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == '_0.001*_z + _1000*_sec - _y'

        new_expr = store.convert_expression_recursively(expr, microsecond)
        assert str(new_expr) == '-_1000*_y + _1e+06*_sec + _z'

    def test_abs_ceil_floor(self, store, sec, ms):
        expr = sp.Abs(sec)
        new_expr = store.convert_expression_recursively(expr, None)
        assert new_expr is expr

        expr = sp.Abs(sec)
        new_expr = store.convert_expression_recursively(expr, ms)
        assert isinstance(new_expr, sp.Abs)
        assert str(new_expr) == 'Abs(_1000*_sec)'
        assert new_expr.args[0].args[1] is sec

        expr = sp.floor(sec)
        new_expr = store.convert_expression_recursively(expr, store.get_unit('second'))
        assert new_expr is expr

        expr = sp.floor(sec)
        new_expr = store.convert_expression_recursively(expr, ms)
        assert isinstance(new_expr, sp.floor)
        assert str(new_expr) == 'floor(_1000*_sec)'
        assert new_expr.args[0].args[1] is sec

        expr = sp.ceiling(sec)
        new_expr = store.convert_expression_recursively(expr, ms)
        assert isinstance(new_expr, sp.ceiling)
        assert str(new_expr) == 'ceiling(_1000*_sec)'
        assert new_expr.args[0].args[1] is sec

    def test_exp_log_trig(self, store, dim, sec, ms):
        dimensionless = store.get_unit('dimensionless')
        z = Variable('z', ms)

        expr = sp.exp(dim)
        assert store.convert_expression_recursively(expr, dimensionless) is expr

        expr = sp.log(sec / z)
        new_expr = store.convert_expression_recursively(expr, None)
        assert isinstance(new_expr, sp.log)
        assert str(new_expr) == 'log(_1000*_sec/_z)'

        expr = sp.sin(z / sec)
        new_expr = store.convert_expression_recursively(expr, dimensionless)
        assert isinstance(new_expr, sp.sin)
        assert str(new_expr) == 'sin(_0.001*_z/_sec)'

    def test_relations(self, store, sec, ms):
        msec = Variable('msec', ms)
        z = Variable('z', ms)

        expr = sec < z
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == '_sec < _0.001*_z'

        expr = z > sec
        new_expr = store.convert_expression_recursively(expr, store.get_unit('dimensionless'))
        assert str(new_expr) == '_z > _1000*_sec'

        # Case with no conversion
        expr = msec < z
        new_expr = store.convert_expression_recursively(expr, None)
        assert new_expr is expr

    def test_piecewise(self, store, ms, microsecond, sec, cm):
        # This also checks more complex nested expressions
        dimensionless = store.get_unit('dimensionless')
        _2 = Quantity(2, dimensionless)
        dim = Variable('dim', dimensionless)
        z = Variable('z', ms)

        expr = sp.Piecewise(
            (sec, sp.And(sec < z, dim > _2)),
            (z, dim < _2),
            (z * _2, True),
        )

        # Units of result will be chosen from first case, i.e. second
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == (
            'Piecewise((_sec, (_2 < _dim) & (_sec < _0.001*_z)), (_0.001*_z, _2 > _dim), (_0.001*_2*_z, True))')

        # Units of result are specified
        new_expr = store.convert_expression_recursively(expr, ms)
        assert str(new_expr) == (
            'Piecewise((_1000*_sec, (_2 < _dim) & (_sec < _0.001*_z)), (_z, _2 > _dim), (_2*_z, True))')

        # A simpler case with no conversion
        _1 = Quantity('1', ms)
        expr = sp.Piecewise((z, dim < _2), (_1, True))
        assert store.convert_expression_recursively(expr, ms) is expr

    def test_assignment(self, store, sec, ms):
        _10 = Quantity(10, ms)
        expr = sp.Eq(sec, _10)
        new_expr = store.convert_expression_recursively(expr, None)
        assert str(new_expr) == 'Eq(_sec, _0.001*_10)'

    def test_evaluate_units_and_fix(self, store, ms):
        s = store.get_unit('second')

        # Note that sympy re-orders the expression into a consistent canonical form, in this case based on the value
        a = Quantity(10, ms)
        b = Quantity(0.1, s)
        expr = a + b
        assert str(expr) == '_0.1 + _10'
        units, new_expr = store.evaluate_units_and_fix(expr)
        assert str(new_expr) == '_0.001*_10 + _0.1'
        assert units is s

        b = Quantity(20, s)
        expr = a + b
        assert str(expr) == '_10 + _20'
        units, new_expr = store.evaluate_units_and_fix(expr)
        assert str(new_expr) == '_10 + _1000*_20'
        assert units is ms

    def test_evaluate_units_and_fix_with_variable(self, store, ms):
        a = Quantity(10, ms)
        b = Variable('b', store.get_unit('second'))
        expr = a + b
        assert str(expr) == '_10 + _b'
        units, new_expr = store.evaluate_units_and_fix(expr)
        assert str(new_expr) == '_10 + _1000*_b'
        assert units is ms

    # Methods below check error cases
    def test_symbol_wrong_dimensions(self, store, mtr):
        with pytest.raises(UnitConversionError, match='from meter to second'):
            store.convert_expression_recursively(mtr, store.get_unit('second'))

    def test_derivative_wrong_dimensions(self, store, mtr):
        t = Variable('t', store.get_unit('second'))
        expr = sp.Derivative(mtr, t)
        with pytest.raises(UnitConversionError, match='Context: trying to convert'):
            store.convert_expression_recursively(expr, store.get_unit('metre'))

    def test_complex_derivative(self, store, mtr, sec):
        t2 = Variable('t2', store.get_unit('second'))
        expr = sp.Derivative(mtr + mtr, sec, t2)
        with pytest.raises(UnexpectedMathUnitsError,
                           match='only support first order derivatives of single variables'):
            store.convert_expression_recursively(expr, store.get_unit('metre'))

    def test_mul_wrong_dimensions(self, store, mtr):
        _1 = Quantity(1, store.get_unit('second'))
        expr = mtr * _1
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, store.get_unit('metre'))

    def test_error_matrix(self, store):
        expr = sp.Matrix([[1, 0], [0, 1]])
        with pytest.raises(UnexpectedMathUnitsError):
            store.convert_expression_recursively(expr, None)

    def test_error_exponent_not_dimensionless(self, store, mtr):
        _1 = Quantity(1, store.get_unit('second'))
        expr = mtr ** _1
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, None)

    def test_error_exponent_not_number(self, store, dim):
        _1 = Quantity(1, store.get_unit('second'))
        expr = _1 ** dim
        with pytest.raises(InputArgumentMustBeNumberError):
            store.convert_expression_recursively(expr, None)

    def test_pow_wrong_dimensions(self, store, sec):
        _2 = Quantity(2, store.get_unit('dimensionless'))
        expr = sec ** _2
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, store.get_unit('second'))

    def test_add_wrong_dimensions(self, store, sec):
        _2 = Quantity(2, store.get_unit('dimensionless'))
        expr = sec + _2
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, None)

    def test_relational_must_be_dimensionless(self, store, sec):
        y = Variable('y', store.get_unit('second'))
        expr = sec < y
        with pytest.raises(BooleanUnitsError):
            store.convert_expression_recursively(expr, store.get_unit('second'))

    def test_relational_dimension_mismatch(self, store, sec):
        y = Variable('y', store.get_unit('metre'))
        expr = sec <= y
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, store.get_unit('dimensionless'))

    def test_piecewise_condition_not_dimensionless(self, store, sec, mtr):
        expr = sp.Piecewise((mtr, sec), (-mtr, True))
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, store.get_unit('metre'))

    def test_piecewise_cannot_convert_result(self, store, mtr, dim):
        expr = sp.Piecewise((mtr, dim), (-mtr, True))
        with pytest.raises(UnitConversionError):
            store.convert_expression_recursively(expr, store.get_unit('second'))

    def test_exp_must_be_dimensionless(self, store, dim):
        expr = sp.exp(dim)
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            store.convert_expression_recursively(expr, store.get_unit('second'))

    def test_number_must_be_dimensionless(self, store):
        with pytest.raises(InputArgumentsMustBeDimensionlessError):
            store.convert_expression_recursively(sp.E, store.get_unit('second'))

    def test_error_unsupported_math(self, store, sec):
        expr = sp.Integral(sec**2, sec)
        with pytest.raises(UnexpectedMathUnitsError):
            store.convert_expression_recursively(expr, None)

    def test_evaluate_units_and_fix_with_unfixable_expr(self, store):
        _10 = Quantity(10, store.get_unit('metre'))
        _20 = Quantity(20, store.get_unit('second'))
        expr = _10 + _20
        with pytest.raises(UnitConversionError, match='Context: trying to evaluate'):
            store.evaluate_units_and_fix(expr)
