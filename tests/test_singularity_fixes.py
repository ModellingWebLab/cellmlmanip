import os

import pytest
from sympy import (
    Abs,
    Float,
    Function,
    I,
    Piecewise,
    exp,
    pi,
)

from cellmlmanip._singularity_fixes import (
    _fix_expr_parts,
    _float_dummies,
    _generate_piecewise,
    _get_singularity,
    _remove_singularities,
    _is_negative_power,
    remove_fixable_singularities,
)
from cellmlmanip.model import Quantity, Variable

from . import shared


OXMETA = 'https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#'


def to_quant(val):
    return Quantity(val, 'dimensionless')


V = Variable(name='V', units='millivolt')
CAL = Variable(name='CAL', units='millivolt')
EXPR1 = 0.4 * V - 18.0
EXPR2 = (-0.16 * V - 1.6) / (exp(-0.16 * V - 1.6) - 1.0)
EXPR3 = -2.1 * (0.06 * V + 1)**2.0
EXPR4 = (-0.26 * V - 2.6) / (exp(-0.26 * V - 2.6) - 1.0)


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


def test_generate_piecewise():
    sp, Vmin, Vmax = 10.0, 9.0, 11.0
    f_Vmin = EXPR1.xreplace({V: Vmin})
    f_Vmax = EXPR1.xreplace({V: Vmax})
    pw = _generate_piecewise(EXPR1, V, sp, Vmin, Vmax)
    assert pw == Piecewise((f_Vmin + ((V - Vmin) / (Vmax - Vmin)) * (f_Vmax - f_Vmin),
                           Abs(V - sp) < Abs((Vmax - Vmin) / 2)), (EXPR1, True))


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
    assert str(sing1) == str(sing2) == '(_-10.0000006250000, _-9.99999937500000, _-10.0000000000000)'


def test_get_singularity2():
    assert str(_get_singularity(EXPR3 / (exp(EXPR3) -1.0), V, 1e-7, exp)) == '(_-16.6666667666667, _-16.6666665666667, _-16.6666666666667)'


def test_get_singularity_no_sing_1():
    assert _get_singularity(EXPR2, V, 1e-7, exp_) == (None, None, None)


def test_get_singularity_no_sing_2():
    assert _get_singularity(CAL / (exp(CAL) -1), V, 1e-7, exp) == (None, None, None)


def test_get_singularity_no_sing_3():
    # test power with other variable still in
    assert _get_singularity(5.0 * V**CAL, V, 1e-7, exp) == (None, None, None)


def test_fix_expr_parts():
    str(_fix_expr_parts(EXPR2, V, 1e-7, exp)) ==\
        ('(_-10.0000006250000, _-9.99999937500000, _-10.0000000000000, '
         '(-0.16*_V - 1.6)/(-1.0 + 0.201896517994655*exp(-0.16*_V)), True)')


def test_fix_expr_parts_add_sing():
    #test adding 2 singularities with same singularity point
    expr = EXPR2 + EXPR4
    str(_fix_expr_parts(expr, V, 1e-7, exp)) ==\
        ('(_-10.0000006250000, _-9.99999937500000, _-10.0000000000000, '
         '(-0.16*_V - 1.6)/(-1.0 + 0.201896517994655*exp(-0.16*_V)) + '
         '(-0.26*_V - 2.6)/(-1.0 + 0.0742735782143339*exp(-0.26*_V)), True)')


def test_singularity_with_variable_in():
    expr =  (V + CAL) / (exp(V + CAL) -1)
    assert str(_get_singularity(expr, V, 1e-7, exp)) == '(_1.00000000000000e-7 - _CAL, _-1.00000000000000e-7 - _CAL, -_CAL)'

def test_multiple_singularities():
    assert str(_get_singularity((EXPR3 + 1.0) / (exp(EXPR3 +1.0) -1.0), V, 1e-7, exp)) == '(_-16.6666667666667, _-16.6666665666667, _-16.6666666666667)'


def test_multiply_singularities():
    expr = EXPR2 * EXPR4
    print(_fix_expr_parts(EXPR2, V, 1e-7, exp) )
    print()
    print(_fix_expr_parts(EXPR4, V, 1e-7, exp) )
    print()
    print(_fix_expr_parts(expr, V, 1e-7, exp) )
    print()
    print(_fix_expr_parts(EXPR4 * EXPR2, V, 1e-7, exp) )
    assert False
    assert str(_remove_singularities(EXPR2, V)) ==\
        ('(True, Piecewise(((-_-10.0000006250000 + _V)*((-0.16*_-9.99999937500000 - 1.6)/(-1.0 + 0.201896517994655*'
         'exp(-0.16*_-9.99999937500000)) - (-0.16*_-10.0000006250000 - 1.6)/(-1.0 + 0.201896517994655*'
         'exp(-0.16*_-10.0000006250000)))/(-_-10.0000006250000 + _-9.99999937500000) + '
         '(-0.16*_-10.0000006250000 - 1.6)/(-1.0 + 0.201896517994655*exp(-0.16*_-10.0000006250000)), '
         'Abs(_-10.0000000000000 - _V) < Abs(_-10.0000006250000/2 - _-9.99999937500000/2)), '
         '((-0.16*_V - 1.6)/(-1.0 + 0.201896517994655*exp(-0.16*_V)), True)))')
    

def test_remove_singularities():
    assert str(_remove_singularities(EXPR2, V)) ==\
        ('(True, Piecewise(((-_-10.0000006250000 + _V)*((-0.16*_-9.99999937500000 - 1.6)/(-1.0 + 0.201896517994655*'
         'exp(-0.16*_-9.99999937500000)) - (-0.16*_-10.0000006250000 - 1.6)/(-1.0 + 0.201896517994655*'
         'exp(-0.16*_-10.0000006250000)))/(-_-10.0000006250000 + _-9.99999937500000) + '
         '(-0.16*_-10.0000006250000 - 1.6)/(-1.0 + 0.201896517994655*exp(-0.16*_-10.0000006250000)), '
         'Abs(_-10.0000000000000 - _V) < Abs(_-10.0000006250000/2 - _-9.99999937500000/2)), '
         '((-0.16*_V - 1.6)/(-1.0 + 0.201896517994655*exp(-0.16*_V)), True)))')


def test_wrong_argument_exclude(model):
    #check piecewises in model without and with singularities fixed
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
