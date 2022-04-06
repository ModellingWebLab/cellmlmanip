import os

import pytest
from sympy import (
    And,
    Float,
    Function,
    Le,
    Or,
    Piecewise,
    Symbol,
    exp,
    log,
)

from cellmlmanip import parser
from cellmlmanip._singularity_fixes import (
    _fix_expr_parts,
    _float_dummies,
    _generate_piecewise,
    _get_singularity,
    _is_negative_power,
    _remove_singularities,
    _solve_real,
    remove_fixable_singularities,
    subs_parsed_math_funcs,
)
from cellmlmanip.model import Quantity, Variable

from . import shared


OXMETA = 'https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#'


def to_quant(val):
    return Quantity(val, 'dimensionless')


V = Variable(name='V', units='millivolt')
CAL = Variable(name='CAL', units='dimensionless')

EXPR1 = 0.4 * V - 18.0
EXPR2 = (-0.16 * V - 1.6) / (exp(-0.16 * V - 1.6) - 1.0)
EXPR3 = -2.1 * (0.06 * V + 1)**2
EXPR4 = (-0.26 * V - 2.6) / (exp(-0.26 * V - 2.6) - 1.0)
EXPR5 = -2.1 * (0.06 * V + 1)**2.0


class exp_(Function):

    def _eval_is_real(self):
        return self.args[0].is_real

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        assert argindex == 1
        return self


@pytest.fixture(scope='session')
def model():
    return shared.load_model('beeler_reuter_model_1977')


def test_subs_parsed_math_funcs():
    v = Symbol('v')
    parser.SIMPLE_MATHML_TO_SYMPY_CLASSES['exp'] = exp_
    expr = 1 / (exp_(v) - 1.0)
    sustituted = subs_parsed_math_funcs(expr)
    parser.SIMPLE_MATHML_TO_SYMPY_CLASSES['exp'] = exp
    assert sustituted == 1 / (exp(v) - 1.0)


def test_generate_piecewise():
    sp, Vmin, Vmax = 10.0, 11.0, 9.0
    pw = _generate_piecewise(EXPR1, V, sp, Vmin, Vmax)
    Vmin, Vmax = Vmax, Vmin  # also checking Vmin/Vmax are sorted correctly

    f_Vmin = EXPR1.xreplace({V: Vmin})
    f_Vmax = EXPR1.xreplace({V: Vmax})
    range_condition = And(Le(Vmin, V), Le(V, Vmax))

    assert pw == Piecewise((f_Vmin + ((V - Vmin) / (Vmax - Vmin)) * (f_Vmax - f_Vmin),
                           range_condition), (EXPR1, True))


def test_generate_piecewise2():
    sp, Vmin, Vmax = 10.0, 9.0 + CAL, 11.0 + CAL
    pw = _generate_piecewise(EXPR1, V, sp, Vmin, Vmax)

    f_Vmin = EXPR1.xreplace({V: Vmin})
    f_Vmax = EXPR1.xreplace({V: Vmax})
    range_condition = Or(And(Le(Vmin, V), Le(V, Vmax)), And(Le(Vmax, V), Le(V, Vmin)))

    assert pw == Piecewise((f_Vmin + ((V - Vmin) / (Vmax - Vmin)) * (f_Vmax - f_Vmin),
                           range_condition), (EXPR1, True))


def test_float_dummies():
    assert str(_float_dummies(12.0 * V)) == str(to_quant(Float(12.0)) * V)


def test_is_negative_power():
    assert not _is_negative_power(EXPR1)
    assert _is_negative_power(1 / EXPR1)
    assert not _is_negative_power(V**2)
    assert _is_negative_power(1 / V**2)
    assert not _is_negative_power(V**CAL)
    assert not _is_negative_power(1 / V**CAL)


def test_get_singularity():
    expr2 = EXPR2.replace(exp, exp_)
    sing1 = _get_singularity(EXPR2, V, 1e-7, exp)
    sing2 = _get_singularity(expr2, V, 1e-7, exp_)
    assert str(sing1) == str(sing2) == '[(_-10.0000006250000, _-9.99999937500000, _-10.0000000000000)]'


def test_get_singularity2():
    assert str(_get_singularity(EXPR3 / (exp(EXPR3) - 1.0), V, 1e-7, exp)) == \
        '[(_-16.6666667666667, _-16.6666665666667, _-16.6666666666667)]'


def test_get_singularity_no_sing_1():
    assert _get_singularity(EXPR2, V, 1e-7, exp_) == []


def test_get_singularity_no_sing_2():
    assert _get_singularity(CAL / (exp(CAL) - 1), V, 1e-7, exp) == []


def test_get_singularity_no_sing_3():
    # test power with other variable still in
    assert _get_singularity(5.0 * V**CAL, V, 1e-7, exp) == []


def test_get_singularity_no_sing_4():
    # test complex expressions with other variables
    # using exp_ to rigger solveset to return an intersection with Reals (it doesn't know the result will be a real)
    expr = (V + log(CAL)) / (exp_(V + log(CAL)) - 1.0)
    assert str(_get_singularity(expr, V, 1e-7, exp_)) == \
        '[(_1.00000000000000e-7 - log(_CAL), _-1.00000000000000e-7 - log(_CAL), -log(_CAL))]'


def test_get_singularity_no_sing_5():
    # test complex expressions with other variables
    # using exp will mean that when matching the top it will look for
    # exp(V + log(CAL)) which will simplify to CAL*exp(V)
    # so it'll find a Z of CAL instead of a number
    expr = (V + log(CAL)) / (exp(V + log(CAL)) - 1.0)
    assert str(_get_singularity(expr, V, 1e-7, exp)) == \
        '[(_1.00000000000000e-7 - log(_CAL), _-1.00000000000000e-7 - log(_CAL), -log(_CAL))]'


def test_singularity_with_variable_in():
    expr = (V + CAL) / (exp(V + CAL) - 1)
    assert str(_get_singularity(expr, V, 1e-7, exp)) == \
        '[(_1.00000000000000e-7 - _CAL, _-1.00000000000000e-7 - _CAL, -_CAL)]'


def test_multiple_singularities_multiple_solutions():
    print(_get_singularity((EXPR3 + 1.0) / (exp(EXPR3 + 1.0) - 1.0), V, 1e-7, exp))
    assert str(_get_singularity((EXPR3 + 1.0) / (exp(EXPR3 + 1.0) - 1.0), V, 1e-7, exp)) == \
        ('[(_-28.1677594223726, _-28.1677592223726, _-28.1677593223726), '
         '(_-5.16557411096076, _-5.16557391096076, _-5.16557401096076)]')


def test_solve_real():
    assert str(set(_solve_real(EXPR3, V))) == '{-16.6666666666667}'
    assert str(set(_solve_real(V + log(CAL), V))) == '{-log(_CAL)}'


def test_multiply_singularities_same_singularity_point_with_var():
    # multiply 2 expressions with the same singularity point, but different Vmin/Vmax
    # due to the variables we can't take the wider range so we get 2 seperate singularities
    sing1 = (V + CAL) / (exp(V + CAL) - 1)
    sing2 = 5 * (V + CAL) / (exp(5 * (V + CAL)) - 1)
    assert str(_get_singularity(sing1 * sing2, V, 1e-7, exp)) == \
        ('[(_2.00000000000000e-8 - _CAL, _-2.00000000000000e-8 - _CAL, -_CAL), '
         '(_1.00000000000000e-7 - _CAL, _-1.00000000000000e-7 - _CAL, -_CAL)]')


def test_multiply_singularities_different_singularity_point():
    sing1 = (V - 0.8) / (exp(V - 0.8) - 1.0)
    assert str(_get_singularity(EXPR2 * sing1, V, 1e-7, exp)) == \
        ('[(_0.800000100000000, _0.799999900000000, _0.800000000000000), '
         '(_-10.0000006250000, _-9.99999937500000, _-10.0000000000000)]')


def test_fix_expr_parts():
    str(_fix_expr_parts(EXPR2, V, 1e-7, exp)) ==\
        ('(_-10.0000006250000, _-9.99999937500000, _-10.0000000000000, '
         '(-0.16*_V - 1.6)/(-1.0 + 0.201896517994655*exp(-0.16*_V)), True)')


def test_fix_expr_parts_add_sing():
    # test adding 2 singularities with same singularity point
    expr = EXPR2 + EXPR4
    str(_fix_expr_parts(expr, V, 1e-7, exp)) ==\
        ('(_-10.0000006250000, _-9.99999937500000, _-10.0000000000000, '
         '(-0.16*_V - 1.6)/(-1.0 + 0.201896517994655*exp(-0.16*_V)) + '
         '(-0.26*_V - 2.6)/(-1.0 + 0.0742735782143339*exp(-0.26*_V)), True)')


def test_multiply_singularities_same_singularity_point():
    assert str(_fix_expr_parts(EXPR2 * EXPR4, V, 1e-7, exp)) ==\
        ('(_-10.0000006250000, _-9.99999937500000, _-10.0000000000000, '
         '(-0.26*_V - 2.6)*(-0.16*_V - 1.6)/((-1.0 + 0.0742735782143339*exp(-0.26*_V))*'
         '(-1.0 + 0.201896517994655*exp(-0.16*_V))), True)')


def test_multiply_singularities_different_singularity_point2():
    sing1 = (V - 0.8) / (exp(V - 0.8) - 1.0)
    expected = open(os.path.join(os.path.dirname(__file__), 'test_singularity.txt'), 'r').read()
    assert str(_fix_expr_parts(EXPR2 * sing1, V, 1e-7, exp)) == expected


def test_remove_singularities():
    expected = open(os.path.join(os.path.dirname(__file__), 'test_singularity2.txt'), 'r').read()
    assert str(_remove_singularities(EXPR2, V)) == expected


def test_wrong_argument_exclude(model):
    # check piecewises in model without and with singularities fixed
    with pytest.raises(TypeError, match='exclude is expected to be a set'):
        model.remove_fixable_singularities(V, exclude=[])


def test_remove_fixable_singularities(model):
    V = model.get_variable_by_ontology_term((OXMETA, 'membrane_voltage'))

    old_piecewises = tuple(str(eq.lhs) for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(old_piecewises) == 1

    remove_fixable_singularities(model, V, [])
    new_piecewises = tuple(eq.lhs for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(new_piecewises) == 3

    additional_piecewises = [v for v in new_piecewises if str(v) not in old_piecewises]
    assert str(additional_piecewises) == '[_sodium_current_m_gate$alpha_m, _time_independent_outward_current$i_K1]'


def test_fix_singularity_equations():
    model = shared.load_model('beeler_reuter_model_1977')
    V = model.get_variable_by_ontology_term((OXMETA, 'membrane_voltage'))

    old_piecewises = tuple(str(eq.lhs) for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(old_piecewises) == 1

    model.remove_fixable_singularities(V)
    new_piecewises = tuple(eq.lhs for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(new_piecewises) == 3

    additional_piecewises = [v for v in new_piecewises if str(v) not in old_piecewises]
    assert str(additional_piecewises) == '[_sodium_current_m_gate$alpha_m, _time_independent_outward_current$i_K1]'
