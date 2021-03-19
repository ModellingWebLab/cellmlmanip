import os

from sympy import (
    Function,
    Piecewise,
    exp,
    log,
)

from cellmlmanip.model import Quantity, Variable
from cellmlmanip.parser import Transpiler
from cellmlmanip.printer import Printer

from . import shared


def to_quant(val):
    return Quantity(val, 'dimensionless')


V = Variable(name='V', units='millivolt')
Z = to_quant(6.2)
U = to_quant(2) * V + to_quant(5.0)
SP = (-(log(Z) + to_quant(5.0)) / to_quant(2))
P = Printer()


class exp_(Function):

    def _eval_is_real(self):
        return self.args[0].is_real

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        assert argindex == 1
        return self


def test_no_V(caplog):
    model = shared.load_model('test_no_V_no_odes')
    model.remove_fixable_singularities()
    assert 'WARNING' in caplog.text
    assert "Has no membrane_voltage tagged, cannot remove fixable singuarities." in caplog.text


def test_fix_singularity_equations():
    # check piecewises in model without and with singularities fixed
    model = shared.load_model('faber_rudy_2000')
    old_piecewises = tuple(eq.lhs for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(old_piecewises) == 14

    model.remove_fixable_singularities()
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
    model = shared.load_model('beeler_reuter_model_1977')

    old_piecewises = tuple(str(eq.lhs) for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(old_piecewises) == 1

    model.remove_fixable_singularities()
    new_piecewises = tuple(eq.lhs for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(new_piecewises) == 3

    additional_piecewises = [v for v in new_piecewises if str(v) not in old_piecewises]
    assert str(additional_piecewises) == '[_sodium_current_m_gate$alpha_m, _time_independent_outward_current$i_K1]'


def test_fix_singularity_equations3():
    # check this still works when the parser is set to generate a differente exp function
    Transpiler.set_mathml_handler('exp', exp_)  # restore exp function
    model = shared.load_model('bondarenko_szigeti_bett_kim_rasmusson_2004_apical')
    old_piecewises = tuple(str(eq.lhs) for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(old_piecewises) == 1
    assert str([eq.rhs for eq in model.equations if str(eq.lhs) ==
                'slow_delayed_rectifier_potassium_current$alpha_n'][0] ==
           '_4.81333e-06*(_26.5 + _membrane$V)/(_1 - exp_(-_0.128*(_26.5 + _membrane$V)))')

    model.remove_fixable_singularities()
    Transpiler.set_mathml_handler('exp', exp)  # restore exp function

    new_piecewises = tuple(eq.lhs for eq in model.equations if eq.rhs.has(Piecewise))
    assert len(new_piecewises) == 2

    additional_piecewises = [v for v in new_piecewises if str(v) not in old_piecewises]
    assert str(additional_piecewises) == '[_slow_delayed_rectifier_potassium_current$alpha_n]'

    expected = open(os.path.join(os.path.dirname(__file__), 'singularity.txt'), 'r').read()

    assert str([eq.rhs for eq in model.equations if str(eq.lhs) ==
                'slow_delayed_rectifier_potassium_current$alpha_n'][0]) == expected
