import os
from xml.dom import pulldom

import sympy

from cellmlmanip import mathml2sympy


class TestParser(object):

    @staticmethod
    def make_mathml(content_xml):
        xml = '<?xml version="1.0"?>' \
              '<math xmlns="http://www.w3.org/1998/Math/MathML">%s</math>' % content_xml
        return xml

    def assert_equal(self, content_xml, sympy_expression):
        mathml_string = self.make_mathml(content_xml)
        transpiled_sympy = mathml2sympy.parse_string(mathml_string)
        # print(mathml_string, "⟶", transpiled_sympy, "==", sympy_expression)
        assert transpiled_sympy == sympy_expression

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

    def test_negative(self):
        self.assert_equal('<apply><minus/><ci>x</ci></apply>',
                          [-sympy.Symbol('x')])

    def test_divide(self):
        self.assert_equal('<apply><divide/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Symbol('a') / sympy.Symbol('b')])

    def test_eq(self):
        self.assert_equal('<apply><eq/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Eq(sympy.Symbol('a'), sympy.Symbol('b'))])

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

    def test_pi(self):
        self.assert_equal('<pi/>',
                          [sympy.pi])

    def test_mod(self):
        self.assert_equal('<apply><rem/><ci>a</ci><ci>b</ci></apply>',
                          [sympy.Symbol('a') % sympy.Symbol('b')])

    def test_log(self):
        self.assert_equal('<apply><log/><ci>x</ci></apply>',
                          [sympy.log(sympy.Symbol('x'), 10)])

    def test_diff(self):
        time = sympy.Symbol('time')
        V = sympy.Function('V')
        self.assert_equal('<apply><diff/><bvar><ci>time</ci></bvar><ci>V</ci></apply>',
                          [sympy.Derivative(V(time), time)])

    def test_piecewise(self):
        x = sympy.Symbol('x')
        self.assert_equal('<piecewise>'
                          '<piece><cn>0</cn><apply><lt/><ci>x</ci><cn>0</cn></apply></piece>'
                          '<otherwise><ci>x</ci></otherwise>'
                          '</piecewise>',
                          [sympy.Piecewise((0, x < 0.0), (x, True))])

    def test_multiple_equalities(self):
        self.assert_equal('<apply><eq/><ci>x</ci><cn>1</cn></apply>'
                          '<apply><eq/><ci>y</ci><cn>2</cn></apply>',
                          [sympy.Eq(sympy.Symbol('x'), 1), sympy.Eq(sympy.Symbol('y'), 2)])

    def test_cellml_namespace(self):
        mathml_xml = '<math xmlns="http://www.w3.org/1998/Math/MathML" ' \
                     'xmlns:cellml="http://www.cellml.org/cellml/1.0#"> <apply><cn ' \
                     'cellml:units="dimensionless">3</cn></apply></math> '
        transpiled_sympy = mathml2sympy.parse_string(mathml_xml)
        assert transpiled_sympy == [sympy.Number(3.0)]

    def test_diff_eq(self):
        t, i_Stim, i_Na, i_K, i_L, Cm = sympy.symbols('time i_Stim i_Na i_K i_L Cm')
        V = sympy.Function('V')
        # From hodgkin_huxley_squid_axon_model_1952_modified.cellml
        self.assert_equal('<apply><eq/><apply><diff/><bvar><ci>time</ci></bvar><ci>V</ci></apply'
                          '><apply><divide/><apply><minus/><apply><plus/><ci>i_Stim</ci><ci>i_Na'
                          '</ci><ci>i_K</ci><ci>i_L</ci></apply></apply><ci>Cm</ci></apply'
                          '></apply>',
                          [sympy.Eq(sympy.Derivative(V(t), t), -(i_Stim+i_Na+i_K+i_L) / Cm)])

    def test_noble_1962(self):
        cellml_path = os.path.join(os.path.dirname(__file__), "noble_model_1962.cellml")

        document = pulldom.parse(cellml_path)
        components = []
        for event, node in document:
            if event == pulldom.START_ELEMENT and node.tagName == 'math':
                document.expandNode(node)
                components.append(node)
        for component in components:
            eqs = mathml2sympy.parse_dom(component)
            for eq in eqs:
                pass
                # print(eq)
