import pytest

from cellmlmanip.units import UnitStore


class TestUnits(object):
    # These represent CellML <units><unit>...</unit></units> elements
    test_definitions = {
        'ms': [{'units': 'second', 'prefix': 'milli'}],
        'per_ms': [{'units': 'ms', 'exponent': '-1'}],
        'usec': [{'units': 'second', 'prefix': 'micro'}],
        'mV': [{'units': 'volt', 'prefix': 'milli'}],
        'per_mV': [{'units': 'volt', 'prefix': 'milli', 'exponent': '-1'}],
        'uV': [{'units': 'volt', 'prefix': 'micro'}],
        'mV_per_ms': [{'units': 'mV', 'exponent': '1'}, {'units': 'ms', 'exponent': '-1'}],
        'mV_per_s': [{'units': 'mV', 'exponent': '1'}, {'units': 'second', 'exponent': '-1'}],
        'mV_per_usec': [{'units': 'mV', 'exponent': '1'},
                        {'prefix': 'micro', 'units': 'second', 'exponent': '-1'}],
        'mM': [{'prefix': 'milli', 'units': 'mole'}, {'units': 'litre', 'exponent': '-1'}],
        'mM_per_ms': [{'units': 'mM'}, {'units': 'ms', 'exponent': '-1'}],
        'milli_mole': [{'prefix': 'milli', 'units': 'mole'}],
        'millisecond': [{'prefix': 'milli', 'units': 'second'}],
        'ms_power_prefix': [{'prefix': '-3', 'units': 'second'}],
    }

    @pytest.fixture(scope="class")
    def quantity_store(self):
        qs = UnitStore(model=None)
        for unit_name, unit_attributes in self.test_definitions.items():
            qs.add_custom_unit(unit_name, unit_attributes)
        return qs

    def test_quantity_translation(self, quantity_store):
        unit_registry = quantity_store.ureg

        assert quantity_store.get_quantity('dimensionless') == unit_registry.dimensionless

        # Units defined in the test CellML <model>:
        for name in TestUnits.test_definitions.keys():
            quantity_store.get_quantity(name)

        # Pint built-in units
        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('kilogram'),
            unit_registry.kilogram
        )

        # Custom units defined in CellML example
        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('per_ms'),
            unit_registry.millisecond**(-1)
        )

        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('mM_per_ms'),
            (unit_registry.milli_mole / unit_registry.liter) / unit_registry.ms
        )

        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('mV_per_usec'),
            unit_registry.mV / unit_registry.usec
        )

        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('ms_power_prefix'),
            unit_registry.millisecond
        )

    def test_conversion_factor(self, quantity_store):
        ureg = quantity_store.ureg
        assert quantity_store.get_conversion_factor(1*ureg.ms, ureg.second) == 0.001
        assert quantity_store.get_conversion_factor(1*ureg.volt, ureg.mV) == 1000.0

        assert quantity_store.get_conversion_factor(
            1 * quantity_store.get_quantity('milli_mole'),
            quantity_store.get_quantity('mole')
        ) == 0.001

    def test_unit_extraction(self, quantity_store):
        pass
        # TODO: test pint unit extraction from sympy expressions
        # kq = (5*units.mile/(2*units.hour + 10*units.minute))**(8*units.gram)
        # assert model.units.summarise_units(eq) == \
        #        (units.mile/(units.hour + units.minute))**units.gram

        # millivolts = units.Quantity('millivolts', units.voltage, units.milli * units.volts, 'mV')
        # x, y = sympy.symbols('x y')
        # eq = (millivolts / units.millisecond)*sympy.Derivative(x, y)
        # assert model.units.summarise_units(eq) == (millivolts / units.milliseconds)
