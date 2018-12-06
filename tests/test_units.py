import pint
import pint.unit
import pytest
import sympy

from cellmlmanip.units import QuantityStore, QuantityStorePints


class TestUnits(object):
    test_definitions = {
        'ms': [{'units': 'second', 'prefix': 'milli'}] ,
        'per_ms': [{'units': 'ms', 'exponent': '-1'}] ,
        'usec': [{'units': 'second', 'prefix': 'micro'}] ,
        'mV': [{'units': 'volt', 'prefix': 'milli'}] ,
        'per_mV': [{'units': 'volt', 'prefix': 'milli', 'exponent': '-1'}] ,
        'uV': [{'units': 'volt', 'prefix': 'micro'}] ,
        'mV_per_ms': [{'units': 'mV', 'exponent': '1'}, {'units': 'ms', 'exponent': '-1'}] ,
        'mV_per_s': [{'units': 'mV', 'exponent': '1'}, {'units': 'second', 'exponent': '-1'}] ,
        'mV_per_usec': [{'units': 'mV', 'exponent': '1'}, {'prefix': 'micro', 'units': 'second', 'exponent': '-1'}] ,
        'mM': [{'prefix': 'milli', 'units': 'mole'}, {'units': 'litre', 'exponent': '-1'}],
        'mM_per_ms': [{'units': 'mM'}, {'units': 'ms', 'exponent': '-1'}],
        'milli_mole': [{'prefix': 'milli', 'units': 'mole'}]
    }

    @pytest.fixture(scope="class")
    def quantity_store(self):
        qs = QuantityStorePints(TestUnits.test_definitions)
        return qs

    @staticmethod
    def check_unit_equivalent(unit1, unit2):
        assert isinstance(unit1, pint.unit._Unit)
        assert isinstance(unit2, pint.unit._Unit)

        assert unit1.dimensionality == unit2.dimensionality

        quantity1 = 1 * unit1
        quantity2 = 1 * unit2
        return quantity1.to(quantity2).magnitude == quantity1.magnitude

    def test_quantity_translation(self, quantity_store):
        unit_registry = quantity_store.ureg

        assert quantity_store.get_quantity('dimensionless') == unit_registry.dimensionless

        # Units defined in the test CellML <model>:
        unit_names = ['ms', 'per_ms', 'usec', 'mV', 'per_mV', 'uV', 'mV_per_ms', 'mV_per_s',
                      'mV_per_usec', 'mM', 'mM_per_ms']
        for name in unit_names:
            quantity_store.get_quantity(name)

        # Pint built-in units
        assert TestUnits.check_unit_equivalent(
            quantity_store.get_quantity('kilogram'),
            unit_registry.kilogram
        )

        # Custom units defined in CellML example
        assert TestUnits.check_unit_equivalent(
            quantity_store.get_quantity('per_ms'),
            unit_registry.millisecond**(-1)
        )

        assert TestUnits.check_unit_equivalent(
            quantity_store.get_quantity('usec'),
            unit_registry.microsecond
        )

        assert TestUnits.check_unit_equivalent(
            quantity_store.get_quantity('mM_per_ms'),
            (unit_registry.millimole / unit_registry.liter) / unit_registry.millisecond
        )

        assert TestUnits.check_unit_equivalent(
            quantity_store.get_quantity('milli_mole'),
            unit_registry.millimole
        )

        assert TestUnits.check_unit_equivalent(
            quantity_store.get_quantity('mV_per_usec'),
            unit_registry.millivolt / unit_registry.microsecond
        )
