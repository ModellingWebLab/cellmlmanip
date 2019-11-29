import pytest


def test_add_preferred_custom_unit_name(simple_ode_model, OXMETA):
    time_var = simple_ode_model.get_symbol_by_ontology_term(OXMETA, "time")
    assert str(simple_ode_model.units.summarise_units(time_var)) == "ms"
    simple_ode_model.units.add_preferred_custom_unit_name('millisecond', [{'prefix': 'milli', 'units': 'second'}])
    assert str(simple_ode_model.units.summarise_units(time_var)) == "millisecond"
    # add_custom_unit does not allow adding already existing units but add_preferred_custom_unit_name does since we
    # cannot know in advance if a model will already have the unit named this way. To test this we add the same unit
    # again
    simple_ode_model.units.add_preferred_custom_unit_name('millisecond', [{'prefix': 'milli', 'units': 'second'}])
    assert str(simple_ode_model.units.summarise_units(time_var)) == "millisecond"


def test_conversion_factor_original(model_simple_units):
    symbol_b1 = model_simple_units.get_symbol_by_cmeta_id("b_1")
    equation = model_simple_units.get_equations_for([symbol_b1])
    factor = model_simple_units.units.get_conversion_factor(
        quantity=1 * model_simple_units.units.summarise_units(equation[0].lhs),
        to_unit=model_simple_units.units.ureg('us').units)
    assert factor == 1000


def test_conversion_factor_bad_types(model_simple_units):
    symbol_b1 = model_simple_units.get_symbol_by_cmeta_id("b_1")
    equation = model_simple_units.get_equations_for([symbol_b1])
    expression = equation[0].lhs
    to_unit = model_simple_units.units.ureg('us').units
    from_unit = model_simple_units.units.summarise_units(expression)
    quantity = 1 * from_unit
    # no source unit
    with pytest.raises(AssertionError, match='^No unit given as source.*'):
        model_simple_units.units.get_conversion_factor(to_unit=to_unit)
    with pytest.raises(AssertionError, match='^No unit given as source.*'):
        model_simple_units.units.get_conversion_factor(to_unit)

    # no target unit
    with pytest.raises(TypeError):
        model_simple_units.units.get_conversion_factor(from_unit=from_unit)
    # multiple sources
    with pytest.raises(AssertionError, match='^Multiple target.*'):
        model_simple_units.units.get_conversion_factor(to_unit, from_unit=from_unit, quantity=quantity)
    # incorrect types
    with pytest.raises(AssertionError, match='^from_unit must be of type pint:Unit$'):
        model_simple_units.units.get_conversion_factor(to_unit, from_unit=quantity)
    with pytest.raises(AssertionError, match='^quantity must be of type pint:Quantity$'):
        model_simple_units.units.get_conversion_factor(to_unit, quantity=from_unit)
    with pytest.raises(AssertionError, match='^expression must be of type Sympy expression$'):
        model_simple_units.units.get_conversion_factor(to_unit, expression=quantity)

    # unit to unit
    assert model_simple_units.units.get_conversion_factor(to_unit=to_unit, from_unit=from_unit) == 1000
    # quantity to unit
    assert model_simple_units.units.get_conversion_factor(to_unit=to_unit, quantity=quantity) == 1000
    # expression to unit
    assert model_simple_units.units.get_conversion_factor(to_unit=to_unit, expression=expression) == 1000


def test_conversion_factor_same_units(model_simple_units):
    symbol_b = model_simple_units.get_symbol_by_cmeta_id("b")
    equation = model_simple_units.get_equations_for([symbol_b])
    expression = equation[1].rhs
    to_unit = model_simple_units.units.ureg('per_ms').units
    from_unit = model_simple_units.units.summarise_units(expression)
    quantity = 1 * from_unit
    # quantity to unit
    assert model_simple_units.units.get_conversion_factor(to_unit=to_unit, quantity=quantity) == 1
    # unit to unit
    assert model_simple_units.units.get_conversion_factor(to_unit=to_unit, from_unit=from_unit) == 1
    # expression to unit
    assert model_simple_units.units.get_conversion_factor(to_unit=to_unit, expression=expression) == 1
