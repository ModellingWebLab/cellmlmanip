import os

import pytest

from cellmlmanip import load_model


OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"


class TestAslanidiModel:
    @pytest.fixture(scope="class")
    def model(self):
        cellml = os.path.join(
            os.path.dirname(__file__),
            "cellml_files",
            "aslanidi_model_2009.cellml"
        )
        return load_model(cellml)

    def test_derivative_symbols(self, model):
        derivs = model.get_derivative_symbols()
        derivs_order = model.get_derivative_symbols(order_by_order_added=True)
        assert len(derivs) == len(derivs_order) == 29

        sorted_derivs = sorted(derivs, key=lambda deriv: str(deriv))
        sorted_derivs_order = sorted(derivs_order, key=lambda deriv: str(deriv))
        for i in range(len(sorted_derivs)):
            assert sorted_derivs[i] == sorted_derivs_order[i]
        print(str(derivs_order))

        assert str(derivs_order) == '[Derivative(_membrane$V, _environment$time), '\
            'Derivative(_sodium_current_m_gate$m, _environment$time), '\
            'Derivative(_sodium_current_h1_gate$h1, _environment$time), '\
            'Derivative(_sodium_current_h2_gate$h2, _environment$time), '\
            'Derivative(_L_type_Ca_channel_d_L_gate$d_L, _environment$time), '\
            'Derivative(_L_type_Ca_channel_f_L_gate$f_L, _environment$time), '\
            'Derivative(_T_type_Ca_channel_d_T_gate$d_T, _environment$time), '\
            'Derivative(_T_type_Ca_channel_f_T_gate$f_T, _environment$time), '\
            'Derivative(_Ca_independent_transient_outward_K_current_r_gate$r, _environment$time), '\
            'Derivative(_Ca_independent_transient_outward_K_current_s1_gate$s1, _environment$time), '\
            'Derivative(_Ca_independent_transient_outward_K_current_s2_gate$s2, _environment$time), '\
            'Derivative(_Ca_independent_transient_outward_K_current_s3_gate$s3, _environment$time), '\
            'Derivative(_delayed_rectifier_K_current_z_gate$z, _environment$time), '\
            'Derivative(_delayed_rectifier_K_current_pa_gate$p_a, _environment$time), '\
            'Derivative(_delayed_rectifier_K_current_pi_gate$p_i, _environment$time), '\
            'Derivative(_intracellular_ion_concentrations$Na_i, _environment$time), '\
            'Derivative(_intracellular_ion_concentrations$Ca_i, _environment$time), '\
            'Derivative(_intracellular_ion_concentrations$K_i, _environment$time), '\
            'Derivative(_intracellular_Ca_buffering$O_C, _environment$time), '\
            'Derivative(_intracellular_Ca_buffering$O_TC, _environment$time), '\
            'Derivative(_intracellular_Ca_buffering$O_TMgC, _environment$time), '\
            'Derivative(_intracellular_Ca_buffering$O_TMgMg, _environment$time), '\
            'Derivative(_cleft_space_ion_concentrations$K_c, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$Ca_rel, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$Ca_up, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$O_Calse, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$F1, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$F2, _environment$time), '\
            'Derivative(_Ca_handling_by_the_SR$F3, _environment$time)]'
