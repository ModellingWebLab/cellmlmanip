import pytest
from sympy import (
    Function,
    Piecewise,
    exp,
    log,
)

from cellmlmanip.model import Quantity, Variable
from cellmlmanip.parser import XmlNs
from cellmlmanip.singularity_fixes import fix_singularity_equations, new_expr

from . import shared


def to_quant(val):
    return Quantity(val, 'dimensionless')


V = Variable(name='V', units='millivolt')
Z = to_quant(6.2)
U = to_quant(2) * V + to_quant(5.0)
SP = (-(log(Z) + to_quant(5.0)) / to_quant(2))


@pytest.fixture(scope='session')
def expr1():
    return (Z * V - Z * SP) / (exp(U) * -Z + to_quant(1.0))


@pytest.fixture(scope='session')
def expr2():
    return (Z * V - Z * SP) / (exp(U) * Z - to_quant(1.0))


class exp_(Function):

    def _eval_is_real(self):
        return self.args[0].is_real

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        assert argindex == 1
        return self


def test_no_singularity():
    expr = (exp_(U) * -Z + to_quant(1.0)) / to_quant(25.0)
    changed, new_ex = new_expr(expr, V)
    assert str(new_expr(expr, V)) == str((False, expr))


def test_no_singularity2():
    expr = to_quant(1) / (to_quant(1) + exp_(to_quant(1.5669291338582676) - to_quant(0.07874015748031496) * V))
    assert str(new_expr(expr, V)) == str((False, expr))


def test_new_expr1(expr1):
    print(str(new_expr(expr1, V)))
    assert str(new_expr(expr1, V)) == ('(True, Piecewise(((-_-3.4122746960255230 + _V)*(-(_-3.4122746960255230*_6.2 - '
                                       '_6.2*(-_5 - log(_6.2))/_2)/(_1 - _6.2*exp(_-3.4122746960255230*_2 + _5)) + (_'
                                       '-3.4122745960255229*_6.2 - _6.2*(-_5 - log(_6.2))/_2)/(_1 - _6.2*exp(_-3.4122'
                                       '745960255229*_2 + _5)))/(_-3.4122745960255229 - _-3.4122746960255230) + (_-3.4'
                                       '122746960255230*_6.2 - _6.2*(-_5 - log(_6.2))/_2)/(_1 - _6.2*exp(_-3.41227469'
                                       '60255230*_2 + _5)), Abs(_-3.4122746460255229 - _V) < Abs(_-3.4122745960255229/'
                                       '2 - _-3.4122746960255230/2)), ((_6.2*_V - _6.2*(-_5 - log(_6.2))/_2)/(_1 - _6.'
                                       '2*exp(_2*_V + _5)), True)))')
    assert str(new_expr(expr1, V, U_offset=1e-5)) == ('(True, Piecewise(((-_-3.4122796460255230 + _V)*(-(_-3.412279646'
                                                      '0255230*_6.2 - _6.2*(-_5 - log(_6.2))/_2)/(_1 - _6.2*exp(_-3.4'
                                                      '122796460255230*_2 + _5)) + (_-3.4122696460255229*_6.2 - _6.2*('
                                                      '-_5 - log(_6.2))/_2)/(_1 - _6.2*exp(_-3.4122696460255229*_2 + '
                                                      '_5)))/(_-3.4122696460255229 - _-3.4122796460255230) + (_-3.4122'
                                                      '796460255230*_6.2 - _6.2*(-_5 - log(_6.2))/_2)/(_1 - _6.2*exp('
                                                      '_-3.4122796460255230*_2 + _5)), Abs(_-3.4122746460255229 - _V) '
                                                      '< Abs(_-3.4122696460255229/2 - _-3.4122796460255230/2)), ((_6.2'
                                                      '*_V - _6.2*(-_5 - log(_6.2))/_2)/(_1 - _6.2*exp(_2*_V + _5)), '
                                                      'True)))')
    assert str(new_expr(expr1, V, exp_function=exp_)) == ('(False, (_6.2*_V - _6.2*(-_5 - log(_6.2))/_2)/(_1 - _6.2*exp'
                                                          '(_2*_V + _5)))')
    expr = expr1.replace(exp, exp_)
    assert str(new_expr(expr, V, exp_function=exp_)) == ('(True, Piecewise(((-_-3.4122746960255230 + _V)*(-(_-3.412274'
                                                         '6960255230*_6.2 - _6.2*(-_5 - log(_6.2))/_2)/(_1 - _6.2*exp_'
                                                         '(_-3.4122746960255230*_2 + _5)) + (_-3.4122745960255229*_6.2'
                                                         ' - _6.2*(-_5 - log(_6.2))/_2)/(_1 - _6.2*exp_(_-3.4122745960'
                                                         '255229*_2 + _5)))/(_-3.4122745960255229 - _-3.41227469602552'
                                                         '30) + (_-3.4122746960255230*_6.2 - _6.2*(-_5 - log(_6.2))/_2'
                                                         ')/(_1 - _6.2*exp_(_-3.4122746960255230*_2 + _5)), Abs(_-3.41'
                                                         '22746460255229 - _V) < Abs(_-3.4122745960255229/2 - _-3.4122'
                                                         '746960255230/2)), ((_6.2*_V - _6.2*(-_5 - log(_6.2))/_2)/(_1'
                                                         ' - _6.2*exp_(_2*_V + _5)), True)))')


def test_new_expr2(expr2):
    assert str(new_expr(expr2, V)) == ('(True, Piecewise(((-_-3.4122746960255230 + _V)*(-(_-3.4122746960255230*_6.2 - '
                                       '_6.2*(-_5 - log(_6.2))/_2)/(-_1 + _6.2*exp(_-3.4122746960255230*_2 + _5)) + ('
                                       '_-3.4122745960255229*_6.2 - _6.2*(-_5 - log(_6.2))/_2)/(-_1 + _6.2*exp(_-3.41'
                                       '22745960255229*_2 + _5)))/(_-3.4122745960255229 - _-3.4122746960255230) + (_-3'
                                       '.4122746960255230*_6.2 - _6.2*(-_5 - log(_6.2))/_2)/(-_1 + _6.2*exp(_-3.41227'
                                       '46960255230*_2 + _5)), Abs(_-3.4122746460255229 - _V) < Abs(_-3.41227459602552'
                                       '29/2 - _-3.4122746960255230/2)), ((_6.2*_V - _6.2*(-_5 - log(_6.2))/_2)/(-_1 +'
                                       ' _6.2*exp(_2*_V + _5)), True)))')


def test_new_expr3(expr1):
    assert str(new_expr(1 / expr1, V)) == ('(True, Piecewise(((-_-3.4122746960255230 + _V)*((_1 - _6.2*exp(_-3.412274'
                                           '5960255229*_2 + _5))/(_-3.4122745960255229*_6.2 - _6.2*(-_5 - log(_6.2))/_'
                                           '2) - (_1 - _6.2*exp(_-3.4122746960255230*_2 + _5))/(_-3.4122746960255230*'
                                           '_6.2 - _6.2*(-_5 - log(_6.2))/_2))/(_-3.4122745960255229 - _-3.41227469602'
                                           '55230) + (_1 - _6.2*exp(_-3.4122746960255230*_2 + _5))/(_-3.4122746960255'
                                           '230*_6.2 - _6.2*(-_5 - log(_6.2))/_2), Abs(_-3.4122746460255229 - _V) < Ab'
                                           's(_-3.4122745960255229/2 - _-3.4122746960255230/2)), ((_1 - _6.2*exp(_2*_'
                                           'V + _5))/(_6.2*_V - _6.2*(-_5 - log(_6.2))/_2), True)))')


def test_new_expr4(expr2):
    assert str(new_expr(1 / expr2, V)) == ('(True, Piecewise(((-_-3.4122746960255230 + _V)*((-_1 + _6.2*exp(_-3.41227'
                                           '45960255229*_2 + _5))/(_-3.4122745960255229*_6.2 - _6.2*(-_5 - log(_6.2))/'
                                           '_2) - (-_1 + _6.2*exp(_-3.4122746960255230*_2 + _5))/(_-3.412274696025523'
                                           '0*_6.2 - _6.2*(-_5 - log(_6.2))/_2))/(_-3.4122745960255229 - _-3.412274696'
                                           '0255230) + (-_1 + _6.2*exp(_-3.4122746960255230*_2 + _5))/(_-3.4122746960'
                                           '255230*_6.2 - _6.2*(-_5 - log(_6.2))/_2), Abs(_-3.4122746460255229 - _V) <'
                                           ' Abs(_-3.4122745960255229/2 - _-3.4122746960255230/2)), ((-_1 + _6.2*exp('
                                           '_2*_V + _5))/(_6.2*_V - _6.2*(-_5 - log(_6.2))/_2), True)))')


def test_new_expr5(expr1):  # Try with numbers in place of Quantity Dummies
    expr = expr1.xreplace({v: float(str(v)) for v in expr1.atoms(Quantity)})
    print(new_expr(expr, V))
    assert str(new_expr(expr, V)) == ('(True, Piecewise(((-_-3.41227469602552 + _V)*(-(6.2*_-3.41227469602552 + 21.156'
                                      '1028053582)/(1.0 - 920.161586435975*exp(2.0*_-3.41227469602552)) + (6.2*_-3.412'
                                      '27459602552 + 21.1561028053582)/(1.0 - 920.161586435975*exp(2.0*_-3.41227459602'
                                      '552)))/(_-3.41227459602552 - _-3.41227469602552) + (6.2*_-3.41227469602552 + 21'
                                      '.1561028053582)/(1.0 - 920.161586435975*exp(2.0*_-3.41227469602552)), Abs(_-3.4'
                                      '1227464602552 - _V) < Abs(_-3.41227459602552/2 - _-3.41227469602552/2)), ((6.2*'
                                      '_V + 21.1561028053582)/(1.0 - 920.161586435975*exp(2.0*_V)), True)))')


def test_fix_singularity_equations():
    # check piecewises in model without and with singularities fixed
    model = shared.load_model('faber_rudy_2000', skip_singularity_fixes=True)
    old_piecewises = tuple(eq.lhs for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(old_piecewises) == 14

    V = model.get_variable_by_ontology_term((XmlNs.OXMETA.value, 'membrane_voltage'))
    fix_singularity_equations(model, V, [])
    new_piecewises = tuple(eq.lhs for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(new_piecewises) == 24

    additional_piecewises = [v for v in new_piecewises if v not in old_piecewises]
    assert str(additional_piecewises) == ('[_L_type_Ca_channel$I_CaCa, '
                                          '_L_type_Ca_channel_d_gate$tau_d, '
                                          '_fast_sodium_current_m_gate$alpha_m, '
                                          '_L_type_Ca_channel$I_CaK, '
                                          '_L_type_Ca_channel$I_CaNa, '
                                          '_non_specific_calcium_activated_current$I_ns_K, '
                                          '_non_specific_calcium_activated_current$I_ns_Na, '
                                          '_rapid_delayed_rectifier_potassium_current_xr_gate$tau_xr, '
                                          '_slow_delayed_rectifier_potassium_current_xs1_gate$tau_xs1, '
                                          '_slow_delayed_rectifier_potassium_current_xs2_gate$tau_xs2]')


def test_fix_singularity_equations2():
    # check fixes via load_model
    model = shared.load_model('bondarenko_szigeti_bett_kim_rasmusson_2004_apical', skip_singularity_fixes=True)
    old_piecewises = tuple(str(eq.lhs) for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(old_piecewises) == 1

    fixed_model = shared.load_model('bondarenko_szigeti_bett_kim_rasmusson_2004_apical', skip_singularity_fixes=False)
    new_piecewises = tuple(eq.lhs for eq in fixed_model.equations if eq.rhs.has(Piecewise))
    assert len(new_piecewises) == 2

    additional_piecewises = [v for v in new_piecewises if str(v) not in old_piecewises]
    assert str(additional_piecewises) == '[_slow_delayed_rectifier_potassium_current$alpha_n]'
