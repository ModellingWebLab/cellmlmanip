"""
Test translating MathML -> SymPy -> callback suitable for SciPy odeint() -> solving with odeint
"""
import collections

import numpy as np
import sympy
from numpy.testing import assert_almost_equal
from scipy.integrate import odeint
from sympy.utilities.lambdify import lambdastr

from cellmlmanip import mathml2sympy


class TestOdeSolving(object):
    """
    Solving Jonathan's examples from https://github.com/ModellingWebLab/cellmlmanip/issues/2
    """
    TIME_POINTS = 100
    TEST_PRECISION = 6  # abs(desired-actual) < 1.5 * 10**(-decimal)

    @staticmethod
    def make_mathml(content_xml):
        """
        Wraps the Content Markup in <math> tags
        """
        xml = '<?xml version="1.0"?>' \
              '<math xmlns="http://www.w3.org/1998/Math/MathML">%s</math>' % content_xml
        return xml

    @staticmethod
    def get_sympy(mathml):
        sympy_expression = mathml2sympy.parse_string(mathml)
        print("SymPy expression(s):", sympy_expression)
        return sympy_expression

    @staticmethod
    def check_results(desired, actual):
        for exact, numerical in zip(desired, actual):
            assert_almost_equal(exact, numerical, decimal=TestOdeSolving.TEST_PRECISION)

    @staticmethod
    def solve(ode_func, initial_conditions, time_points):
        solution, infodict = odeint(ode_func, initial_conditions, time_points, full_output=True)
        print('odeint() message:', infodict['message'])
        return solution

    @staticmethod
    def dummify_undefined_functions(expr):
        """SymPy Derivative dummification code from https://stackoverflow.com/a/29940130"""
        mapping = {}

        # Replace all Derivative terms
        for der in expr.atoms(sympy.Derivative):
            atom_name = der.expr.func.__name__
            var_names = [var.name for var in der.variables]
            name = "d%s_d%s" % (atom_name, 'd'.join(var_names))
            mapping[der] = sympy.Symbol(name)

        # Replace undefined functions
        from sympy.core.function import AppliedUndef
        for atom in expr.atoms(AppliedUndef):
            atom_name = atom.func.__name__
            mapping[atom] = sympy.Symbol(atom_name)

        return expr.subs(mapping)

    @staticmethod
    def sympy_to_callback(sympy_expressions):
        """
        Takes a list of SymPy expressions that specifies a system of ODEs. Each equality *must*
        follow the form:

            sympy.Eq(sympy.Derivative(sympy.Function, <SymPy expression>))

        This method will extract a list of symbols and generate a SciPy odeint()-compliant
        function. The order of symbols in the list is the order that must be used when passing
        initial conditions in odeint()

        :param sympy_expressions: list of SymPy expressions
        :return: lambda_symbols, odeint_wrapper: tuple of ([sympy.Symbol], Callback)
        """

        derivative_lambdas = list()

        for expr in sympy_expressions:
            if isinstance(expr, sympy.Eq) and isinstance(expr.lhs, sympy.Derivative):
                # Keep an ordered set of symbols used in this derivative
                expr_symbols = collections.OrderedDict()

                # First symbol is always the symbol for the function of this derivative
                derivative_symbol = TestOdeSolving.dummify_undefined_functions(expr.lhs.args[0])
                expr_symbols[derivative_symbol] = None

                # Add any other symbols used on the RHS of the equation
                for symbol in expr.rhs.free_symbols:
                    expr_symbols[symbol] = None

                # Convert to list
                expr_symbols = list(expr_symbols.keys())

                # Create lambda function for this expression
                print('%s -> %s' % (expr, lambdastr(expr_symbols, expr.rhs)))
                expr_lambda = sympy.lambdify(expr_symbols, expr.rhs)
                derivative_lambdas.append((expr_symbols, expr_lambda))

        # Get function symbols for which we have derivatives (first tuple, first item in list)
        function_symbols = [x[0][0] for x in derivative_lambdas]
        print('Functions:', function_symbols)

        def odeint_wrapper(symbols_in, _):
            """Callback to be used by the SciPy odeint() function"""
            # map odeint list input to mapping of { symbol: input }
            symbol_in_mapping = dict()
            for index, function_symbol in enumerate(function_symbols):
                symbol_in_mapping[function_symbol] = symbols_in[index]
            # evaluate each derivative lambda, passing in relevant arguments
            lambda_results = list()
            for derivative_symbols, derivative_lambda in derivative_lambdas:
                lambda_args = [symbol_in_mapping[x] for x in derivative_symbols]
                lambda_result = derivative_lambda(*lambda_args)
                lambda_results.append(lambda_result)
            # return list of evaluated derivatives
            return lambda_results

        return function_symbols, odeint_wrapper

    @staticmethod
    def extract_solutions(variables, solution):
        """
        Returns a dictionary mapping SymPy symbols (of functions) to the results column from odeint
        """
        extracted = dict()
        for index, var in enumerate(variables):
            extracted[var] = solution[:, index]
        return extracted

    def test_example_1(self):
        """
        1. d
           ──(V(t)) = 1.0
           dt

        V(0) = 1.0

        Solution: V(t) = t + 1
        """

        mathml = self.make_mathml("""
        <apply>
           <eq />
           <apply>
              <diff />
              <bvar>
                 <ci>t</ci>
              </bvar>
              <ci>V</ci>
           </apply>
           <apply>
              <cn>1</cn>
           </apply>
        </apply>""")

        sympy_expression = self.get_sympy(mathml)

        function_order, odeint_func = self.sympy_to_callback(sympy_expression)

        # Assume we have initial conditions from somewhere else (e.g. cellml file)
        initial_map = {sympy.Symbol('V'): 1.0}
        ode_initial = [x.subs(initial_map) for x in function_order]

        time_points = np.linspace(0, 5., self.TIME_POINTS)

        solution = self.solve(odeint_func, ode_initial, time_points)
        points = self.extract_solutions(function_order, solution)

        self.check_results([t + 1 for t in time_points], points[sympy.Symbol('V')])

    def test_example_2(self):
        """
        2. ⎡ d              d            ⎤
           ⎢ ──(x(t)) = 1,  ──(y(t)) = -1⎥
           ⎣ dt             dt           ⎦

        x(0) = y(0) = 1

        Solution: x(t) = t + 1,  y(t) = 1 - t
        """
        mathml = self.make_mathml("""
        <apply>
           <eq />
           <apply>
              <diff />
              <bvar>
                 <ci>t</ci>
              </bvar>
              <ci>x</ci>
           </apply>
           <apply>
              <cn>1</cn>
           </apply>
        </apply>
        <apply>
           <eq />
           <apply>
              <diff />
              <bvar>
                 <ci>t</ci>
              </bvar>
              <ci>y</ci>
           </apply>
           <apply>
              <cn>-1</cn>
           </apply>
        </apply>""")

        sympy_expression = self.get_sympy(mathml)

        function_order, odeint_func = self.sympy_to_callback(sympy_expression)

        initial_map = {sympy.Symbol('x'): 1.0, sympy.Symbol('y'): 1.0}
        ode_initial = [x.subs(initial_map) for x in function_order]

        time_points = np.linspace(0, 5., self.TIME_POINTS)

        solution = self.solve(odeint_func, ode_initial, time_points)
        points = self.extract_solutions(function_order, solution)

        self.check_results([t + 1 for t in time_points], points[sympy.Symbol('x')])
        self.check_results([1 - t for t in time_points], points[sympy.Symbol('y')])

    def test_example_3(self):
        """
        3. ⎡ d               d           ⎤
           ⎢ ──(x(t)) = -y,  ──(y(t)) = x⎥
           ⎣ dt              dt          ⎦

        x(0) = 0
        y(0) = 1

        Solution: x(t) = -sin(t), y(t) = cos(t)
        """

        mathml = self.make_mathml("""
        <apply>
           <eq />
           <apply>
              <diff />
              <bvar>
                 <ci>t</ci>
              </bvar>
              <ci>x</ci>
           </apply>
           <apply>
              <minus />
              <ci>y</ci>
           </apply>
        </apply>
        <apply>
           <eq />
           <apply>
              <diff />
              <bvar>
                 <ci>t</ci>
              </bvar>
              <ci>y</ci>
           </apply>
           <apply>
              <ci>x</ci>
           </apply>
        </apply>""")

        sympy_expression = self.get_sympy(mathml)

        lambda_symbols, odeint_func = self.sympy_to_callback(sympy_expression)

        # Imagine we have picked up the initial conditions from somewhere else
        initial_map = {sympy.Symbol('x'): 0, sympy.Symbol('y'): 1}
        ode_initial = [x.subs(initial_map) for x in lambda_symbols]

        # Time points for a circle (0, 2 * pi)
        time_points = np.linspace(0, float(2 * sympy.pi), self.TIME_POINTS)

        solution = self.solve(odeint_func, ode_initial, time_points)
        points = self.extract_solutions(lambda_symbols, solution)

        self.check_results([-np.sin(t) for t in time_points], points[sympy.Symbol('x')])
        self.check_results([np.cos(t) for t in time_points], points[sympy.Symbol('y')])

    def test_example_4(self):
        """
        4. ⎡           d                d            ⎤
           ⎢ t₁ = t₂, ───(x(t₁)) = -y, ───(y(t₂)) = x⎥
           ⎣          dt₁              dt₂           ⎦

        """
        pass


if __name__ == '__main__':
    solving = TestOdeSolving()
    print(solving.test_example_1.__doc__)
    solving.test_example_1()
    print(solving.test_example_2.__doc__)
    solving.test_example_2()
    print(solving.test_example_3.__doc__)
    solving.test_example_3()
