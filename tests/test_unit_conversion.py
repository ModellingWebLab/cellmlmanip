import os

import pytest

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
    factor = simple_model.units.get_conversion_factor(quantity=1 * simple_model.units.summarise_units(equation[0].lhs),
                                                      to_unit=simple_model.units.ureg('us').units)
    assert factor == 1000


def test_conversion_factor_bad_types(simple_model):
    simple_model.get_equation_graph(True)  # set up the graph - it is not automatic
    symbol_b1 = simple_model.get_symbol_by_cmeta_id("b_1")
    equation = simple_model.get_equations_for([symbol_b1])
    expression = equation[0].lhs
    to_unit = simple_model.units.ureg('us').units
    from_unit = simple_model.units.summarise_units(expression)
    quantity = 1 * from_unit
    # no source unit
    with pytest.raises(AssertionError, match='^No unit given as source.*'):
        simple_model.units.get_conversion_factor(to_unit=to_unit)
    with pytest.raises(AssertionError, match='^No unit given as source.*'):
        simple_model.units.get_conversion_factor(to_unit)

    # no target unit
    with pytest.raises(TypeError):
        simple_model.units.get_conversion_factor(from_unit=from_unit)
    # multiple sources
    with pytest.raises(AssertionError, match='^Multiple target.*'):
        simple_model.units.get_conversion_factor(to_unit, from_unit=from_unit, quantity=quantity)
    # incorrect types
    with pytest.raises(AssertionError, match='^from_unit must be of type pint:Unit$'):
        simple_model.units.get_conversion_factor(to_unit, from_unit=quantity)
    with pytest.raises(AssertionError, match='^quantity must be of type pint:Quantity$'):
        simple_model.units.get_conversion_factor(to_unit, quantity=from_unit)
    with pytest.raises(AssertionError, match='^expression must be of type Sympy expression$'):
        simple_model.units.get_conversion_factor(to_unit, expression=quantity)

    # unit to unit
    assert simple_model.units.get_conversion_factor(to_unit=to_unit, from_unit=from_unit) == 1000
    # quantity to unit
    assert simple_model.units.get_conversion_factor(to_unit=to_unit, quantity=quantity) == 1000
    # expression to unit
    assert simple_model.units.get_conversion_factor(to_unit=to_unit, expression=expression) == 1000


def test_conversion_factor_same_units(simple_model):
    simple_model.get_equation_graph(True)  # set up the graph - it is not automatic
    symbol_b = simple_model.get_symbol_by_cmeta_id("b")
    equation = simple_model.get_equations_for([symbol_b])
    expression = equation[1].rhs
    to_unit = simple_model.units.ureg('per_ms').units
    from_unit = simple_model.units.summarise_units(expression)
    quantity = 1 * from_unit
    # quantity to unit
    assert simple_model.units.get_conversion_factor(to_unit=to_unit, quantity=quantity) == 1
    # unit to unit
    assert simple_model.units.get_conversion_factor(to_unit=to_unit, from_unit=from_unit) == 1
    # expression to unit
    assert simple_model.units.get_conversion_factor(to_unit=to_unit, expression=expression) == 1

