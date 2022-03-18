import pytest
import sympy

from cellmlmanip.parser import Transpiler


class TestTranspiler(object):
    """
    Tests conversion from mathml to sympy using the Transpiler class.
    """

    @staticmethod
    def make_mathml(content_xml):
        xml = '<?xml version="1.0"?>' \
              '<math xmlns="http://www.w3.org/1998/Math/MathML" ' \
              'xmlns:cellml="http://www.cellml.org/cellml/1.0#">%s</math>' % content_xml
        return xml

    def assert_equal(self, content_xml, sympy_expression):
        transpiler = Transpiler()
        mathml_string = self.make_mathml(content_xml)
        transpiled_sympy = transpiler.parse_string(mathml_string)
        assert transpiled_sympy == sympy_expression

    def assert_raises(self, content_xml, exception_class, match=None):
        transpiler = Transpiler()
        mathml_string = self.make_mathml(content_xml)
        with pytest.raises(exception_class, match=match):
            transpiler.parse_string(mathml_string)

    def test_symbol(self):
        self.assert_equal('<ci>x</ci>',
                          [sympy.Symbol('x')])

    def test_number(self):
        self.assert_equal('<cn>1</cn>',
                          [sympy.Number(1)])

    def test_ignore_comment(self):
        self.assert_equal('<!-- ignore this -->', [])

    def test_ignore_processing(self):
        self.assert_equal('<?PITarget PIContent?>', [])

    def test_plus(self):
        self.assert_equal('<apply><plus/><ci>x</ci><ci>y</ci><ci>z</ci></apply>',
                          [sympy.Symbol('x') + sympy.Symbol('y') + sympy.Symbol('z')])

    def test_mul(self):
        self.assert_equal('<apply><times/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Symbol('a') * sympy.Symbol('b')])

    def test_minus(self):
        self.assert_equal('<apply><minus/><ci>x</ci><ci>y</ci></apply>',
                          [sympy.Symbol('x') - sympy.Symbol('y')])

    def test_minus_nary(self):
        with pytest.raises(TypeError):
            self.assert_equal('<apply><minus/><ci>x</ci><ci>y</ci><ci>z</ci></apply>',
                              [sympy.Symbol('x') - sympy.Symbol('y') - sympy.Symbol('z')])

    def test_negative(self):
        self.assert_equal('<apply><minus/><ci>x</ci></apply>',
                          [-sympy.Symbol('x')])

    def test_divide(self):
        self.assert_equal('<apply><divide/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Symbol('a') / sympy.Symbol('b')])

    def test_divide_nary(self):
        with pytest.raises(TypeError):
            self.assert_equal('<apply><divide/><ci>a</ci><ci>b</ci><ci>c</ci></apply>',
                              [sympy.Symbol('a') / sympy.Symbol('b') / sympy.Symbol('c')])

    def test_eq(self):
        self.assert_equal('<apply><eq/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Eq(sympy.Symbol('a'), sympy.Symbol('b'))])

    def test_neq(self):
        self.assert_equal('<apply><neq/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Ne(sympy.Symbol('a'), sympy.Symbol('b'))])

    def test_leq(self):
        self.assert_equal('<apply><leq/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Symbol('a') <= sympy.Symbol('b')])

    def test_lt(self):
        self.assert_equal('<apply><lt/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Symbol('a') < sympy.Symbol('b')])

    def test_floor(self):
        self.assert_equal('<apply><floor/><ci>a</ci></apply>',
                          [sympy.floor(sympy.Symbol('a'))])

    def test_geq(self):
        self.assert_equal('<apply><geq/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Symbol('a') >= sympy.Symbol('b')])

    def test_gt(self):
        self.assert_equal('<apply><gt/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Symbol('a') > sympy.Symbol('b')])

    def test_and(self):
        self.assert_equal('<apply><and/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Symbol('a') & sympy.Symbol('b')])

    def test_or(self):
        self.assert_equal('<apply><or/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Symbol('a') | sympy.Symbol('b')])

    def test_exp(self):
        self.assert_equal('<apply><exp/><ci>x</ci></apply>',
                          [sympy.exp(sympy.Symbol('x'))])

    def test_power(self):
        self.assert_equal('<apply><power/><ci>x</ci><cn>3</cn></apply>',
                          [sympy.Symbol('x') ** 3.0])

    def test_ln(self):
        self.assert_equal('<apply><ln/><ci>a</ci></apply>',
                          [sympy.ln(sympy.Symbol('a'))])

    def test_abs(self):
        self.assert_equal('<apply><abs/><ci>x</ci></apply>',
                          [sympy.Abs(sympy.Symbol('x'))])

    def test_root(self):
        self.assert_equal('<apply><root/><ci>a</ci></apply>',
                          [sympy.sqrt(sympy.Symbol('a'))])

    def test_root_degree(self):
        self.assert_equal('<apply><root/><degree><ci>n</ci></degree><ci>a</ci></apply>',
                          [sympy.root(sympy.Symbol('a'), sympy.Symbol('n'))])

    def test_degree_error(self):
        self.assert_raises('<degree></degree>', ValueError, match='Expected single value')

    def test_pi(self):
        self.assert_equal('<pi/>',
                          [sympy.pi])

    def test_mod(self):
        self.assert_equal('<apply><rem/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Symbol('a') % sympy.Symbol('b')])

    def test_log(self):
        self.assert_equal('<apply><log/><ci>x</ci></apply>',
                          [sympy.log(sympy.Symbol('x'), 10)])

    def test_log_with_base(self):
        # numeric base
        self.assert_equal('<apply><log/><logbase><cn>3</cn></logbase><ci>x</ci></apply>',
                          [sympy.log(sympy.Symbol('x'), sympy.Float(3.0))])
        # symbolic base
        self.assert_equal('<apply><log/><logbase><ci>y</ci></logbase><ci>x</ci></apply>',
                          [sympy.log(sympy.Symbol('x'), sympy.Symbol('y'))])

    def test_diff(self):
        time = sympy.Symbol('time')
        V = sympy.Symbol('V')
        self.assert_equal('<apply><diff/><bvar><ci>time</ci></bvar><ci>V</ci></apply>',
                          [sympy.Derivative(V, time)])

    def test_diff_with_order(self):
        time = sympy.Symbol('time')
        V = sympy.Symbol('V')
        self.assert_equal('<apply>'
                          '<diff/>'
                          '<bvar><ci>time</ci><degree><cn>2</cn></degree></bvar>'
                          '<ci>V</ci>'
                          '</apply>',
                          [sympy.Derivative(V, time, 2)])

    def test_bvar_error(self):
        self.assert_raises('<bvar/>', ValueError, match='Do not know how to handle')

    def test_piecewise(self):
        x = sympy.Symbol('x')
        self.assert_equal('<piecewise>'
                          '<piece><cn>0</cn><apply><lt/><ci>x</ci><cn>0</cn></apply></piece>'
                          '<otherwise><ci>x</ci></otherwise>'
                          '</piecewise>',
                          [sympy.Piecewise((0.0, x < 0), (x, True))])

        self.assert_equal('<piecewise>'
                          '<piece><cn>10</cn><apply><gt/><ci>x</ci><cn>0</cn></apply></piece>'
                          '<piece><cn>20</cn><apply><gt/><ci>x</ci><cn>1</cn></apply></piece>'
                          '<piece><cn>30</cn><apply><gt/><ci>x</ci><cn>2</cn></apply></piece>'
                          '<otherwise><cn>0</cn></otherwise>'
                          '</piecewise>',
                          [sympy.Piecewise((10.0, x > 0.), (20.0, x > 1.0), (30.0, x > 2.0), (0.0, True))])

    def test_piece_error(self):
        self.assert_raises('<piece><cn>0</cn></piece>', ValueError, match='Need exactly 2 children')

    def test_otherwise_error(self):
        self.assert_raises('<otherwise><cn>0</cn><cn>1</cn></otherwise>', ValueError, match='More than 1 child')

    def test_multiple_relations(self):
        from sympy.abc import (
            a,
            b,
            c,
            x,
            y,
            z,
        )
        eq_xml = '<apply><eq/><ci>x</ci><ci>y</ci><ci>z</ci><cn>2.0</cn></apply>'
        lt_xml = '<apply><lt/><ci>a</ci><ci>b</ci><ci>c</ci><cn>2.0</cn></apply>'
        ge_xml = '<apply><geq/><ci>x</ci><ci>y</ci><ci>z</ci><cn>2.0</cn></apply>'

        self.assert_equal(eq_xml, [sympy.And(sympy.Eq(x, y), sympy.Eq(y, z), sympy.Eq(z, 2.0))])
        self.assert_equal(lt_xml, [sympy.And(sympy.Lt(a, b), sympy.Lt(b, c), sympy.Lt(c, 2.0))])
        self.assert_equal(ge_xml, [sympy.And(sympy.Ge(x, y), sympy.Ge(y, z), sympy.Ge(z, 2.0))])
        self.assert_equal('<apply><and/>%s%s</apply>' % (eq_xml, lt_xml),
                          [sympy.And(
                              sympy.And(sympy.Eq(x, y), sympy.Eq(y, z), sympy.Eq(z, 2.0)),
                              sympy.And(sympy.Lt(a, b), sympy.Lt(b, c), sympy.Lt(c, 2.0))
                          )])

    def test_cellml_namespace(self):
        mathml_xml = '<math xmlns="http://www.w3.org/1998/Math/MathML" ' \
                     'xmlns:cellml="http://www.cellml.org/cellml/1.0#"> <apply><cn ' \
                     'cellml:units="dimensionless">3</cn></apply></math> '
        transpiled_sympy = Transpiler().parse_string(mathml_xml)
        assert transpiled_sympy == [sympy.Number(3.0)]

    def test_diff_eq(self):
        t, i_Stim, i_Na, i_K, i_L, Cm = sympy.symbols('time i_Stim i_Na i_K i_L Cm')
        V = sympy.Symbol('V')
        # From hodgkin_huxley_squid_axon_model_1952_modified.cellml
        self.assert_equal('<apply><eq/><apply><diff/><bvar><ci>time</ci></bvar><ci>V</ci></apply'
                          '><apply><divide/><apply><minus/><apply><plus/><ci>i_Stim</ci><ci>i_Na'
                          '</ci><ci>i_K</ci><ci>i_L</ci></apply></apply><ci>Cm</ci></apply'
                          '></apply>',
                          [sympy.Eq(sympy.Derivative(V, t), -(i_Stim + i_Na + i_K + i_L) / Cm)])

    def test_scientific_notation(self):
        self.assert_equal('<cn type="e-notation">1.234<sep/>5</cn>', [sympy.Number(1.234e5)])

    def test_bad_cn_element(self):
        self.assert_raises('<cn type="e-notation">1.23<sep/>5<sep/>6</cn>', ValueError, match='significand')
        self.assert_raises('<cn type="complex-cartesian"/>', ValueError, match='Unimplemented type attribute')

    def test_xor(self):
        self.assert_equal('<apply><xor/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Xor(sympy.Symbol('a'), sympy.Symbol('b'))])

    def test_not(self):
        self.assert_equal('<apply><not/><ci>a</ci></apply>',
                          [sympy.Not(sympy.Symbol('a'))])

    def test_ceiling(self):
        self.assert_equal('<apply><ceiling/><ci>a</ci></apply>',
                          [sympy.ceiling(sympy.Symbol('a'))])

    def test_min(self):
        self.assert_equal('<apply><min/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Min(sympy.Symbol('a'), sympy.Symbol('b'))])

    def test_max(self):
        self.assert_equal('<apply><max/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Max(sympy.Symbol('a'), sympy.Symbol('b'))])

    def test_trig(self):
        self.assert_equal('<apply><sin/><ci>x</ci></apply>', [sympy.sin(sympy.Symbol('x'))])
        self.assert_equal('<apply><cos/><ci>x</ci></apply>', [sympy.cos(sympy.Symbol('x'))])
        self.assert_equal('<apply><tan/><ci>x</ci></apply>', [sympy.tan(sympy.Symbol('x'))])

        self.assert_equal('<apply><arcsin/><ci>x</ci></apply>', [sympy.asin(sympy.Symbol('x'))])
        self.assert_equal('<apply><arccos/><ci>x</ci></apply>', [sympy.acos(sympy.Symbol('x'))])
        self.assert_equal('<apply><arctan/><ci>x</ci></apply>', [sympy.atan(sympy.Symbol('x'))])

        self.assert_equal('<apply><sinh/><ci>x</ci></apply>', [sympy.sinh(sympy.Symbol('x'))])
        self.assert_equal('<apply><cosh/><ci>x</ci></apply>', [sympy.cosh(sympy.Symbol('x'))])
        self.assert_equal('<apply><tanh/><ci>x</ci></apply>', [sympy.tanh(sympy.Symbol('x'))])

        self.assert_equal('<apply><arcsinh/><ci>x</ci></apply>', [sympy.asinh(sympy.Symbol('x'))])
        self.assert_equal('<apply><arccosh/><ci>x</ci></apply>', [sympy.acosh(sympy.Symbol('x'))])
        self.assert_equal('<apply><arctanh/><ci>x</ci></apply>', [sympy.atanh(sympy.Symbol('x'))])

        self.assert_equal('<apply><sec/><ci>x</ci></apply>', [sympy.sec(sympy.Symbol('x'))])
        self.assert_equal('<apply><csc/><ci>x</ci></apply>', [sympy.csc(sympy.Symbol('x'))])
        self.assert_equal('<apply><cot/><ci>x</ci></apply>', [sympy.cot(sympy.Symbol('x'))])

        self.assert_equal('<apply><arcsec/><ci>x</ci></apply>', [sympy.asec(sympy.Symbol('x'))])
        self.assert_equal('<apply><arccsc/><ci>x</ci></apply>', [sympy.acsc(sympy.Symbol('x'))])
        self.assert_equal('<apply><arccot/><ci>x</ci></apply>', [sympy.acot(sympy.Symbol('x'))])

        self.assert_equal('<apply><sech/><ci>x</ci></apply>', [sympy.sech(sympy.Symbol('x'))])
        self.assert_equal('<apply><csch/><ci>x</ci></apply>', [sympy.csch(sympy.Symbol('x'))])
        self.assert_equal('<apply><coth/><ci>x</ci></apply>', [sympy.coth(sympy.Symbol('x'))])

        self.assert_equal('<apply><arcsech/><ci>x</ci></apply>', [sympy.asech(sympy.Symbol('x'))])
        self.assert_equal('<apply><arccsch/><ci>x</ci></apply>', [sympy.acsch(sympy.Symbol('x'))])
        self.assert_equal('<apply><arccoth/><ci>x</ci></apply>', [sympy.acoth(sympy.Symbol('x'))])

    def test_cn_units(self):
        """ Test if the transpiler uses the number_generator to make number objects. """

        mathml = self.make_mathml('<apply>'
                                  '<eq/>'
                                  '<cn cellml:units="s">2.0</cn>'
                                  '<cn cellml:units="ms">2000.0</cn>'
                                  '</apply>')
        numbers = []

        def number_generator(number, units):
            numbers.append((number, units))
            return sympy.Rational(number)

        transpiler = Transpiler(number_generator=number_generator)
        expr = transpiler.parse_string(mathml)
        for number in expr[0].free_symbols:
            assert isinstance(number, sympy.Rational)

        assert len(numbers) == 2
        assert numbers[0][0] == 2
        assert numbers[1][0] == 2000

    def test_change_handler(self):
        class _exp(sympy.Function):
            pass
        Transpiler.set_mathml_handler('exp', _exp)
        self.assert_equal('<apply><exp/><ci>x</ci></apply>',
                          [_exp(sympy.Symbol('x'))])

    def test_unhandled_tag(self):
        self.assert_raises('<matrix/>', ValueError, match='No handler for element')
