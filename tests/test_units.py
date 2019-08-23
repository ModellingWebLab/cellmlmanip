import pytest
import sympy as sp

from cellmlmanip.model import MetaDummy
from cellmlmanip.units import (ExpressionWithUnitPrinter, UnitCalculator,
                               UnitStore)


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
        'ms_with_multiplier': [{'multiplier': 0.001, 'units': 'second'}],
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
            # what is being tested ??

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

        assert quantity_store.is_unit_equal(
            quantity_store.get_quantity('ms_with_multiplier'),
            unit_registry.millisecond
        )

    def test_conversion_factor(self, quantity_store):
        ureg = quantity_store.ureg
        assert quantity_store.get_conversion_factor(1 * ureg.ms, ureg.second) == 0.001
        assert quantity_store.get_conversion_factor(1 * ureg.volt, ureg.mV) == 1000.0

        assert quantity_store.get_conversion_factor(
            1 * quantity_store.get_quantity('milli_mole'),
            quantity_store.get_quantity('mole')
        ) == 0.001

    def test_add_custom_unit(self):
        unitstore = UnitStore(model=None)
        assert (unitstore._is_unit_defined('newunit') is False)
        # add a new unit called 'newunit' which consists of
        # (2 * second) - I know not realistic
        unit_attributes = [{'multiplier': 2, 'units': 'second'}]
        unitstore.add_custom_unit('newunit', unit_attributes)
        assert (unitstore._is_unit_defined('newunit') is True)

    def test_add_custom_unit_1(self):
        unitstore = UnitStore(model=None)
        assert (unitstore._is_unit_defined('newunit') is False)
        # add a new unit called 'newunit' which consists of
        # (millimetres)
        unit_attributes = [{'prefix': 'milli', 'units': 'metre'}]
        unitstore.add_custom_unit('newunit', unit_attributes)
        assert (unitstore._is_unit_defined('newunit') is True)

    def test_add_custom_unit_2(self):
        unitstore = UnitStore(model=None)
        assert (unitstore._is_unit_defined('newunit') is False)
        # add a new unit called 'newunit' which consists of
        # (milliVolts)
        unit_attributes = [{'prefix': -3, 'units': 'volt'}]
        unitstore.add_custom_unit('newunit', unit_attributes)
        assert (unitstore._is_unit_defined('newunit') is True)

    def test_add_custom_unit_3(self):
        unitstore = UnitStore(model=None)
        assert (unitstore._is_unit_defined('newunit') is False)
        # add a new unit called 'newunit' which consists of
        # (metre per second)
        unit_attributes = [{'units': 'metre'}, {'units': 'second', 'exponent': -1}]
        unitstore.add_custom_unit('newunit', unit_attributes)
        assert (unitstore._is_unit_defined('newunit') is True)

    def test_add_custom_unit_existing(self):
        unitstore = UnitStore(model=None)
        unit_attributes = [{'multiplier': 1, 'units': 'second'}]
        with pytest.raises(AssertionError):
            unitstore.add_custom_unit('second', unit_attributes)
            pytest.fail("Cannot redefine CellML unit <second>")
            pass

    def test_is_unit_defined(self, quantity_store):
        assert (quantity_store._is_unit_defined('ms') is True)
        assert (quantity_store._is_unit_defined('not_a_unit') is False)

    def test_make_pint_unit_definition(self, quantity_store):
        unit_attributes = [{'prefix': -3, 'units': 'metre'},
                           {'prefix': 'milli', 'exponent': -1, 'units': 'second'},
                           {'multiplier': 2, 'units': 'kilograms'}
                           ]
        assert(quantity_store._make_pint_unit_definition('kg_mm_per_ms', unit_attributes) ==
               'kg_mm_per_ms=(meter * 1e-3)*(((second * 0.001))**-1)*(2 * kilogram)')

    # Test UnitCalculator class

    def test_unit_calculator(self, quantity_store):
        ureg = quantity_store.ureg
        a, b, c, d, x, y, z, _1, _2 = [sp.Dummy(x)
                                       for x in ['a', 'b', 'c', 'd', 'x', 'y', 'z', '1', '2']]

        symbol_info = {
            a: MetaDummy('a', ureg.meter, a),
            b: MetaDummy('b', ureg.second, b),
            c: MetaDummy('c', ureg.gram, c),
            d: MetaDummy('d', ureg.meter, d),
            x: MetaDummy('x', ureg.kilogram, x),
            y: MetaDummy('y', ureg.volt, y),
            z: MetaDummy('z', ureg.ampere, z),
            _1: MetaDummy('_1', ureg.kelvin, _1, number=sp.Float(1.0)),
            _2: MetaDummy('_2', ureg.dimensionless, _2, number=sp.Integer(2)),
        }

        unit_calculator = UnitCalculator(ureg, symbol_info)

        assert unit_calculator.traverse(a + a + a + a).units == ureg.meter
        assert unit_calculator.traverse(a + 2*a + 3*a + 4*a).units == ureg.meter  # noqa: E226
        assert unit_calculator.traverse((a * a) / b).units == ureg.meter**2 / ureg.second
        assert unit_calculator.traverse(a**_2).units == ureg.meter**2
        assert unit_calculator.traverse(sp.sqrt(c ** 2)).units == ureg.gram
        assert unit_calculator.traverse(sp.Abs(-2 * y)).units == ureg.volt
        assert unit_calculator.traverse(sp.Piecewise((a, x < 1),
                                                     (a + a, x > 1),
                                                     (3 * a, True))).units == ureg.meter
        assert unit_calculator.traverse(sp.floor(x)).units == ureg.kilogram
        result = unit_calculator.traverse(sp.floor(12.5) * a)
        assert result.units == ureg.meter and result.magnitude == 12.0 * a

        assert unit_calculator.traverse(sp.sqrt(a * d)).units == ureg.meter

        # bad unit expressions
        assert unit_calculator.traverse(sp.exp(3 * c)) is None
        assert unit_calculator.traverse(a + b + c) is None
        assert unit_calculator.traverse(a ** _1) is None
        assert unit_calculator.traverse(sp.Piecewise((a, x < 1), (b, x > 1), (c, True))) is None

    def test_expression_printer(self, quantity_store):
        ureg = quantity_store.ureg
        a, b, c, x, y, z, _1, _2 = [sp.Dummy(x)
                                    for x in ['a', 'b', 'c', 'x', 'y', 'z', '1', '2']]
        symbol_info = {
            a: MetaDummy('a', ureg.meter, a),
            b: MetaDummy('b', ureg.second, b),
            c: MetaDummy('c', ureg.gram, c),
            x: MetaDummy('x', ureg.kilogram, x),
            y: MetaDummy('y', ureg.volt, y),
            z: MetaDummy('z', ureg.ampere, z),
            _1: MetaDummy('_1', ureg.kelvin, _1, number=sp.Float(1.0)),
            _2: MetaDummy('_2', ureg.dimensionless, _2, number=sp.Integer(2)),
        }
        printer = ExpressionWithUnitPrinter(symbol_info)
        assert printer.doprint(a * a) == 'a[meter]**2'
        assert printer.doprint(a / b) == 'a[meter]/b[second]'
        assert printer.doprint(_2 * y) == '2.000000[dimensionless]*y[volt]'
        assert printer.doprint(sp.Derivative(a, b)) == 'Derivative(a[meter], b[second])'

