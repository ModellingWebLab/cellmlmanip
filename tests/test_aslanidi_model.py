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

    def test_initial_value_capacitance(self, model):
        model.get_equation_graph()
        membrane_capacitance = model.get_symbol_by_ontology_term(OXMETA, "membrane_capacitance")
        # was raising KeyError but should not
        assert(model.get_initial_value(membrane_capacitance) == 0.00005)

    def test_initial_value_voltage(self, model):
        model.get_equation_graph()
        membrane_voltage = model.get_symbol_by_ontology_term(OXMETA, "membrane_voltage")
        assert(model.get_initial_value(membrane_voltage) == -80)

    def test_get_derivative_symbols(self, model):
        derivs = model.get_state_symbols()
        assert str(derivs[0]) == '_membrane$V'
        assert str(derivs[1]) == '_sodium_current_m_gate$m'
        assert str(derivs[2]) == '_sodium_current_h1_gate$h1'
        assert str(derivs[3]) == '_sodium_current_h2_gate$h2'
        assert str(derivs[4]) == '_L_type_Ca_channel_d_L_gate$d_L'
        assert str(derivs[5]) == '_L_type_Ca_channel_f_L_gate$f_L'
        assert str(derivs[6]) == '_T_type_Ca_channel_d_T_gate$d_T'
        assert str(derivs[7]) == '_T_type_Ca_channel_f_T_gate$f_T'
        assert str(derivs[8]) == '_Ca_independent_transient_outward_K_current_r_gate$r'
        assert str(derivs[9]) == '_Ca_independent_transient_outward_K_current_s1_gate$s1'
        assert str(derivs[10]) == '_Ca_independent_transient_outward_K_current_s2_gate$s2'
        assert str(derivs[11]) == '_Ca_independent_transient_outward_K_current_s3_gate$s3'
        assert str(derivs[12]) == '_delayed_rectifier_K_current_z_gate$z'
        assert str(derivs[13]) == '_delayed_rectifier_K_current_pa_gate$p_a'
        assert str(derivs[14]) == '_delayed_rectifier_K_current_pi_gate$p_i'
        assert str(derivs[15]) == '_intracellular_ion_concentrations$Na_i'
        assert str(derivs[16]) == '_intracellular_ion_concentrations$K_i'
        assert str(derivs[17]) == '_intracellular_ion_concentrations$Ca_i'
        assert str(derivs[18]) == '_intracellular_Ca_buffering$O_C'
        assert str(derivs[19]) == '_intracellular_Ca_buffering$O_TC'
        assert str(derivs[20]) == '_intracellular_Ca_buffering$O_TMgC'
        assert str(derivs[21]) == '_intracellular_Ca_buffering$O_TMgMg'
        assert str(derivs[22]) == '_cleft_space_ion_concentrations$K_c'
        assert str(derivs[23]) == '_Ca_handling_by_the_SR$O_Calse'
        assert str(derivs[24]) == '_Ca_handling_by_the_SR$Ca_rel'
        assert str(derivs[25]) == '_Ca_handling_by_the_SR$Ca_up'
        assert str(derivs[26]) == '_Ca_handling_by_the_SR$F1'
        assert str(derivs[27]) == '_Ca_handling_by_the_SR$F2'
        assert str(derivs[28]) == '_Ca_handling_by_the_SR$F3'

        derivs = model.get_state_symbols(order_by_order_added=True)
        assert str(derivs[0]) == '_membrane$V'
        assert str(derivs[1]) == '_sodium_current_m_gate$m'
        assert str(derivs[2]) == '_sodium_current_h1_gate$h1'
        assert str(derivs[3]) == '_sodium_current_h2_gate$h2'
        assert str(derivs[4]) == '_L_type_Ca_channel_d_L_gate$d_L'
        assert str(derivs[5]) == '_L_type_Ca_channel_f_L_gate$f_L'
        assert str(derivs[6]) == '_T_type_Ca_channel_d_T_gate$d_T'
        assert str(derivs[7]) == '_T_type_Ca_channel_f_T_gate$f_T'
        assert str(derivs[8]) == '_Ca_independent_transient_outward_K_current_r_gate$r'
        assert str(derivs[9]) == '_Ca_independent_transient_outward_K_current_s1_gate$s1'
        assert str(derivs[10]) == '_Ca_independent_transient_outward_K_current_s2_gate$s2'
        assert str(derivs[11]) == '_Ca_independent_transient_outward_K_current_s3_gate$s3'
        assert str(derivs[12]) == '_delayed_rectifier_K_current_z_gate$z'
        assert str(derivs[13]) == '_delayed_rectifier_K_current_pa_gate$p_a'
        assert str(derivs[14]) == '_delayed_rectifier_K_current_pi_gate$p_i'
        assert str(derivs[15]) == '_intracellular_ion_concentrations$Na_i'
        assert str(derivs[16]) == '_intracellular_ion_concentrations$Ca_i'
        assert str(derivs[17]) == '_intracellular_ion_concentrations$K_i'
        assert str(derivs[18]) == '_intracellular_Ca_buffering$O_C'
        assert str(derivs[19]) == '_intracellular_Ca_buffering$O_TC'
        assert str(derivs[20]) == '_intracellular_Ca_buffering$O_TMgC'
        assert str(derivs[21]) == '_intracellular_Ca_buffering$O_TMgMg'
        assert str(derivs[22]) == '_cleft_space_ion_concentrations$K_c'
        assert str(derivs[23]) == '_Ca_handling_by_the_SR$Ca_rel'
        assert str(derivs[24]) == '_Ca_handling_by_the_SR$Ca_up'
        assert str(derivs[25]) == '_Ca_handling_by_the_SR$O_Calse'
        assert str(derivs[26]) == '_Ca_handling_by_the_SR$F1'
        assert str(derivs[27]) == '_Ca_handling_by_the_SR$F2'
        assert str(derivs[28]) == '_Ca_handling_by_the_SR$F3'
