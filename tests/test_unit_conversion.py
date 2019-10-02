import pytest
import os
import cellmlmanip

OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"


@pytest.fixture
def model():
    return cellmlmanip.load_model(os.path.join(os.path.dirname(__file__), 'cellml_files', "test_simple_odes.cellml"))


@pytest.fixture
def simple_model():
    return cellmlmanip.load_model(os.path.join(os.path.dirname(__file__), 'cellml_files', "simple_model_units.cellml"))


def test_add_preferred_custom_unit_name(model):
    time_var = model.get_symbol_by_ontology_term(OXMETA, "time")
    assert str(model.units.summarise_units(time_var)) == "ms"
    model.units.add_preferred_custom_unit_name('millisecond', [{'prefix': 'milli', 'units': 'second'}])
    assert str(model.units.summarise_units(time_var)) == "millisecond"
    # add_custom_unit does not allow adding already existing units but add_preferred_custom_unit_name does since we
    # cannot know in advance if a model will already have the unit named this way. To test this we add the same unit
    # again
    model.units.add_preferred_custom_unit_name('millisecond', [{'prefix': 'milli', 'units': 'second'}])
    assert str(model.units.summarise_units(time_var)) == "millisecond"


def test_conversion_factor_original(simple_model):
    simple_model.get_equation_graph(True)  # set up the graph - it is not automatic
    symbol_b1 = simple_model.get_symbol_by_cmeta_id("b_1")
    equation = simple_model.get_equations_for([symbol_b1])
    factor = simple_model.units.get_conversion_factor(1 * simple_model.units.summarise_units(equation[0].lhs),
                                                      simple_model.units.ureg('us').units)
    assert factor == 1000

def test_convers_factor(simple_model):
    simple_model.get_equation_graph(True)  # set up the graph - it is not automatic
    symbol_b1 = simple_model.get_symbol_by_cmeta_id("b_1")
    equation = simple_model.get_equations_for([symbol_b1])
    expression = equation[0].lhs
    to_unit = simple_model.units.ureg('us').units
    from_unit = simple_model.units.summarise_units(expression)
    quantity = 1 * from_unit
    assert simple_model.units.get_convers_factor(to_unit=to_unit, quantity=quantity) == 1000
