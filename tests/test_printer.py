#
# Tests conversion of Sympy expressions to Python code
#
import logging
import math

import pytest
import sympy as sp

import cellmlmanip.printer


# Show more logging output
logging.getLogger().setLevel(logging.INFO)


class TestPrinter(object):

    @pytest.fixture(scope="class")
    def x(self):
        return sp.symbols('x')

    @pytest.fixture(scope="class")
    def y(self):
        return sp.symbols('y')

    @pytest.fixture(scope="class")
    def z(self):
        return sp.symbols('z')

    @pytest.fixture(scope="class")
    def p(self):
        return cellmlmanip.printer.Printer()

    def test_numbers(self, p, x):
        # Number types
        assert p.doprint(1) == '1'                  # int
        assert p.doprint(1.2) == '1.2'              # float, short format
        assert p.doprint(math.pi) == '3.141592653589793'  # float, long format
        assert p.doprint(1.436432635636e-123) == '1.436432635636e-123'
        assert p.doprint(x - x) == '0'              # Zero
        assert p.doprint(x / x) == '1'              # One
        assert p.doprint(-x / x) == '-1'            # Negative one
        assert p.doprint(5 * (x / x)) == '5'        # Sympy integer
        assert p.doprint(5.5 * (x / x)) == '5.5'        # Sympy float
        assert p.doprint(sp.Rational(5, 7)) == '5 / 7'  # Sympy rational

        # Special numbers
        assert p.doprint(sp.pi) == 'math.pi'
        assert p.doprint(sp.E) == 'math.e'

    def test_symbols(self, p, x, y):
        # Symbols
        assert p.doprint(x) == 'x'

        # Derivatives
        assert p.doprint(sp.Derivative(x, y)) == 'Derivative(x, y)'

        # Symbol function
        def symbol_function(symbol):
            return symbol.name.upper()

        q = cellmlmanip.printer.Printer(symbol_function)
        assert q.doprint(x) == 'X'
        assert q.doprint(sp.Derivative(x, y)) == p.doprint(sp.Derivative(x, y))

        # Derivative function
        def derivative_function(deriv):
            a = deriv.expr
            b = deriv.variables[0]
            return 'd' + a.name + '/' + 'd' + b.name.upper()

        q = cellmlmanip.printer.Printer(derivative_function=derivative_function)
        assert q.doprint(sp.Derivative(x, y)) == 'dx/dY'
        assert q.doprint(x) == p.doprint(x)

        # Both
        q = cellmlmanip.printer.Printer(symbol_function, derivative_function)
        assert q.doprint(x) == 'X'
        assert q.doprint(sp.Derivative(x, y)) == 'dx/dY'

    def test_addition(self, p, x, y, z):

        # Addition and subtraction
        assert p.doprint(x + y) == 'x + y'
        assert p.doprint(x + y + z) == 'x + y + z'
        assert p.doprint(x - y) == 'x - y'
        assert p.doprint(2 + z) == '2 + z'
        assert p.doprint(z + 2) == '2 + z'
        assert p.doprint(z - 2) == '-2 + z'
        assert p.doprint(2 - z) == '2 - z'
        assert p.doprint(-x) == '-x'
        assert p.doprint(-x - 2) == '-2 - x'

    def test_multiplication(self, p, x, y, z):

        # Multiplication and division
        assert p.doprint(x * y) == 'x * y'
        assert p.doprint(x * y * z) == 'x * y * z'
        assert p.doprint(x / y) == 'x / y'
        assert p.doprint(2 * z) == '2 * z'
        assert p.doprint(z * 5) == '5 * z'
        assert p.doprint(4 / z) == '4 / z'
        assert p.doprint(z / 3) == 'z / 3'
        assert p.doprint(1 / x) == '1 / x'  # Uses pow
        assert p.doprint(1 / (x * y)) == '1 / (x * y)'
        assert p.doprint(1 / -(x * y)) == '-1 / (x * y)'
        assert p.doprint(x + (y + z)) == 'x + y + z'
        assert p.doprint(x * (y + z)) == 'x * (y + z)'
        assert p.doprint(x * y * z) == 'x * y * z'
        assert p.doprint(x + y > x * z), 'x + y > x * z'
        assert p.doprint(x**2 + y**2) == 'x**2 + y**2'
        assert p.doprint(x**2 + 3 * y**2) == 'x**2 + 3 * y**2'
        assert p.doprint(x**(2 + y**2)) == 'x**(2 + y**2)'
        assert p.doprint(x**(2 + 3 * y**2)) == 'x**(2 + 3 * y**2)'
        assert p.doprint(x**-1 * y**-1) == '1 / (x * y)'
        assert p.doprint(x / y / z) == 'x / (y * z)'
        assert p.doprint(x / y * z) == 'x * z / y'
        assert p.doprint(x / (y * z)) == 'x / (y * z)'
        assert p.doprint(x * y**(-2 / (3 * x / x))) == 'x / y**(2 / 3)'

        # Sympy issue #14160
        d = sp.Mul(
            -2,
            x,
            sp.Pow(sp.Mul(y, y, evaluate=False), -1, evaluate=False),
            evaluate=False
        )
        assert p.doprint(d) == '-2 * x / (y * y)'

    def test_powers(self, p, x, y, z):

        # Powers and square roots
        assert p.doprint(sp.sqrt(2)) == 'math.sqrt(2)'
        assert p.doprint(1 / sp.sqrt(2)) == 'math.sqrt(2) / 2'
        assert p.doprint(
            sp.Mul(1, 1 / sp.sqrt(2), eval=False)) == 'math.sqrt(2) / 2'
        assert p.doprint(sp.sqrt(x)) == 'math.sqrt(x)'
        assert p.doprint(1 / sp.sqrt(x)) == '1 / math.sqrt(x)'
        assert p.doprint(x**(x / (2 * x))) == 'math.sqrt(x)'
        assert p.doprint(x**(x / (-2 * x))) == '1 / math.sqrt(x)'
        assert p.doprint(x**-1) == '1 / x'
        assert p.doprint(x**0.5) == 'x**0.5'
        assert p.doprint(x**-0.5) == 'x**(-0.5)'
        assert p.doprint(x**(1 + y)) == 'x**(1 + y)'
        assert p.doprint(x**-(1 + y)) == 'x**(-1 - y)'
        assert p.doprint((x + z)**-(1 + y)) == '(x + z)**(-1 - y)'
        assert p.doprint(x**-2) == 'x**(-2)'
        assert p.doprint(x**3.2) == 'x**3.2'

    def test_trig_functions(self, p, x):

        # Trig functions
        assert p.doprint(sp.acos(x)) == 'math.acos(x)'
        assert p.doprint(sp.acosh(x)) == 'math.acosh(x)'
        assert p.doprint(sp.asin(x)) == 'math.asin(x)'
        assert p.doprint(sp.asinh(x)) == 'math.asinh(x)'
        assert p.doprint(sp.atan(x)) == 'math.atan(x)'
        assert p.doprint(sp.atanh(x)) == 'math.atanh(x)'
        assert p.doprint(sp.ceiling(x)) == 'math.ceil(x)'
        assert p.doprint(sp.cos(x)) == 'math.cos(x)'
        assert p.doprint(sp.cosh(x)) == 'math.cosh(x)'
        assert p.doprint(sp.factorial(x)) == 'math.factorial(x)'
        assert p.doprint(sp.floor(x)) == 'math.floor(x)'
        assert p.doprint(sp.log(x)) == 'math.log(x)'
        assert p.doprint(sp.sin(x)) == 'math.sin(x)'
        assert p.doprint(sp.sinh(x)) == 'math.sinh(x)'
        assert p.doprint(sp.tan(x)) == 'math.tan(x)'
        assert p.doprint(sp.tanh(x)) == 'math.tanh(x)'

    def test_conditions(self, p, x, y, z):

        # Conditions
        assert p.doprint(sp.Eq(x, y)) == 'x == y'
        assert p.doprint(sp.Eq(x, sp.Eq(y, z))) == 'x == (y == z)'
        assert p.doprint(sp.Eq(sp.Eq(x, y), z)) == '(x == y) == z'
        assert p.doprint(sp.Ne(x, y)) == 'x != y'
        assert p.doprint(sp.Gt(x, y)) == 'x > y'
        assert p.doprint(sp.Lt(x, y)) == 'x < y'
        assert p.doprint(sp.Ge(x, y)) == 'x >= y'
        assert p.doprint(sp.Le(x, y)) == 'x <= y'
        e = sp.Eq(sp.Eq(x, 3), sp.Eq(y, 5))
        assert p.doprint(e) == '(x == 3) == (y == 5)'

    def test_boolean_logic(self, p, x, y, z):

        # Boolean logic
        assert p.doprint(True) == 'True'
        assert p.doprint(False) == 'False'
        assert p.doprint(sp.Eq(x, x)) == 'True'
        assert p.doprint(sp.Ne(x, x)) == 'False'
        assert (
            p.doprint(sp.And(sp.Eq(x, y), sp.Eq(x, z))) == 'x == y and x == z')
        assert (
            p.doprint(sp.And(sp.Eq(x, y), sp.Eq(x, z), sp.Eq(x, 2))) ==
            'x == 2 and x == y and x == z')
        assert p.doprint(sp.Or(sp.Eq(x, y), sp.Eq(x, z))) == 'x == y or x == z'
        assert (
            p.doprint(sp.Or(sp.Eq(x, y), sp.Eq(x, z), sp.Eq(x, 2))) ==
            'x == 2 or x == y or x == z')
        a, b, c = x > 2, x > y, x > z
        assert p.doprint(a & b) == 'x > 2 and x > y'
        # 1 or (0 and 0) = 1 = 1 or 0 and 0 -- and binds stronger
        # (1 or 0) and 0 = 0
        assert p.doprint(a | (b & c)) == 'x > 2 or x > y and x > z'
        assert p.doprint((a | b) & c) == 'x > z and (x > 2 or x > y)'

    def test_piecewise_expressions(self, p, x):

        # Piecewise expressions
        e = sp.Piecewise((0, x > 0), (1, x > 1), (2, True))
        assert (
            p.doprint(e) == '((0) if (x > 0) else ((1) if (x > 1) else (2)))')
        e = sp.Piecewise((0, x > 0), (1, x > 1), (2, True), (3, x > 3))
        assert (
            p.doprint(e) == '((0) if (x > 0) else ((1) if (x > 1) else (2)))')
        # Sympy filters out False statements
        e = sp.Piecewise(
            (0, x > 0), (1, x != x), (2, True), (3, x > 3),
            evaluate=False)
        assert p.doprint(e) == '((0) if (x > 0) else (2))'

        # First condition false, multiple true clauses
        e = sp.Piecewise((6, False), (0, x < 5), (1, True), (2, True))
        assert p.doprint(e) == '((0) if (x < 5) else (1))'

        # Middle condition only true
        e = sp.Piecewise((0, x < 5), (1, True), (2, x > 7))
        assert p.doprint(e) == '((0) if (x < 5) else (1))'

    def test_long_expression(self, p, x, y, z):

        # Longer expressions
        assert (
            p.doprint((x + y) / (2 + z / sp.exp(x - y))) ==
            '(x + y) / (2 + z * math.exp(y - x))')
        assert p.doprint((y + sp.sin(x))**-1) == '1 / (y + math.sin(x))'

    def test_unsupported_sympy_items(self, p, x):

        # Unsupported sympy item
        e = sp.Matrix()
        with pytest.raises(ValueError):
            p.doprint(e)

        # Unsupported sympy function
        e = sp.gamma(x)
        with pytest.raises(ValueError):
            p.doprint(e)

    def test_abs(self, p, x, y):
        assert p.doprint(sp.Abs(x + y)) == 'abs(x + y)'
        assert p.doprint(sp.Abs(3.2, evaluate=False)) == 'abs(3.2)'
        assert p.doprint(sp.Abs(-3, evaluate=False)) == 'abs(-3)'

    def test_print_model_equations(model, simple_ode_model, p):
        printed_equations = [p.doprint(eq.lhs) + ' = ' + p.doprint(eq.rhs) for eq in simple_ode_model.equations]
        assert str(printed_equations) == \
            ("['time_units_conversion1$time = 0.001 * environment$time', 'time_units_conversion2$time = 1000.0 * envir"
             "onment$time', 'Derivative(_single_independent_ode$sv1, _environment$time) = 1.0', 'Derivative(_single_od"
             "e_rhs_const_var$sv1, _environment$time) = single_ode_rhs_const_var$a', 'Derivative(_single_ode_rhs_compu"
             "ted_var$sv1, _environment$time) = single_ode_rhs_computed_var$a', 'single_ode_rhs_computed_var$a = -1.0'"
             ", 'derived_from_state_var$dbl_sv1 = 2.0 * single_ode_rhs_computed_var$sv1', 'deriv_on_rhs$sv1_rate = Der"
             "ivative(_single_ode_rhs_computed_var$sv1, _environment$time)', 'Derivative(_circle_x_source$x, _environm"
             "ent$time) = -1.0 * circle_y_implementation$y', 'circle_x_sibling$x2 = 2.0 * circle_x_source$x', 'Derivat"
             "ive(_circle_y_implementation$y, _environment$time) = circle_y_implementation$rhs', 'circle_y_implementat"
             "ion$rhs = 1.0 * circle_x_source$x', 'circle_sibling$local_complex_maths = 1.0 * math.exp(circle_y_implem"
             "entation$y / 2.0) + 5.0 * circle_y_implementation$y / 3.0 + circle_y_implementation$y', 'Derivative(_tim"
             "e_units_conversion1$sv1, _environment$time) = 0.001', 'deriv_on_rhs1a$sv1_rate = Derivative(_time_units_"
             "conversion1$sv1, _environment$time)', 'Derivative(_time_units_conversion2$sv1, _environment$time) = 1000"
             ".0', 'deriv_on_rhs2a$sv1_rate = Derivative(_time_units_conversion2$sv1, _environment$time)', 'Derivative"
             "(_state_units_conversion1$sv1, _environment$time) = 0.001', 'deriv_on_rhs1b$sv1_rate = Derivative(_state"
             "_units_conversion1$sv1, _environment$time)', 'Derivative(_state_units_conversion2$sv1, _environment$time"
             ") = 1000.0', 'deriv_on_rhs2b$sv1_rate = Derivative(_state_units_conversion2$sv1, _environment$time)', 's"
             "ingle_ode_rhs_const_var$a = 1.0']")

    def test_pow_non_commutative(self, p, x, y):
        x1 = sp.Symbol('x', commutative=False)
        y1 = sp.Symbol('y', commutative=False)
        expr = sp.Pow(x1 / y1, -1)
        assert p.doprint(expr) == "(x * y**(-1))**(-1)"
        expr = sp.Pow(x / y, -1)
        assert p.doprint(expr) == "y / x"
