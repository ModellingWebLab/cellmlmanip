"""
Unit handling for CellML models, using the Pint unit library (replaces previous
Sympy-units implementation
"""
import logging
import math
import re

import pint
import sympy
from pint.quantity import _Quantity as Quantity
from pint.unit import _Unit as Unit
from sympy.printing.lambdarepr import LambdaPrinter


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# The full list of supported CellML units
# Taken from https://www.cellml.org/specifications/cellml_1.1/#sec_units
# Some are not defined by Sympy, see comments
CELLML_UNITS = {
    # Base SI units
    'ampere',
    'candela',
    'kelvin',
    'kilogram',
    'meter',
    'mole',
    'second',

    # Derived SI units
    'becquerel',
    'celsius',
    'coulomb',
    'farad',
    'gray',
    'henry',
    'hertz',
    'joule',
    'katal',  # see __add_custom_units()
    'lumen',
    'lux',
    'newton',
    'ohm',
    'pascal',
    'radian',
    'siemens',
    'sievert',
    'steradian',
    'tesla',
    'volt',
    'watt',
    'weber',

    # Convenience units
    'dimensionless',
    'gram',
    'liter',

    # Aliases
    'metre',
    'litre',
}


class UnitStore(object):
    """
    Wraps the underlying Pint UnitRegistry to provide unit handling for the model's Sympy
    expressions. Getting and checking units is handled by this class.
    """
    def __init__(self, model, cellml_def=None):
        """Initialise a QuantityStore instance
        :param cellml_def: a dictionary of <units> definitions from the CellML model. See parser
        for format, essentially: {'name of unit': { [ unit attributes ], [ unit attributes ] } }
        """
        # Initialise the unit registry
        # TODO: create Pint unit definition file for CellML
        self.ureg: pint.UnitRegistry = pint.UnitRegistry()

        # Add default CellML units not provided by Pint
        self.__add_undefined_units()

        # Hold on to custom unit definitions
        self.cellml_definitions = cellml_def if cellml_def else {}
        logger.debug('Found %d CellML unit definitions', len(self.cellml_definitions))

        # CellML units that we've defined in unit registry because of call to get_quantity()
        self.cellml_defined = set()

        # Keep reference to the underlying model, to look up the 'dummy_info' dictionary
        self.model = model

        # Prints Sympy expression in terms of units
        self.printer = PintUnitPrinter(self)

    def __add_undefined_units(self):
        """Adds units required by CellML but not provided by Pint."""
        self.ureg.define('katal = mol / second = kat')

    def get_quantity(self, unit_name: str):
        """Given the name of the unit, this will either (i) return a Quantity from the internal
        store as its already been resolved previously (ii) return the Quantity from Pint if it a
        built-in name or (iii) construct and return a new Quantity object using the
        <units><unit></units> definition in the CellML <model>
        """
        # If this unit is a custom CellML definition and we haven't defined it
        if unit_name in self.cellml_definitions and unit_name not in self.cellml_defined:
            # TODO create cellml pint unit definitions file
            try:
                getattr(self.ureg, unit_name)
            except pint.UndefinedUnitError:
                # Create the unit definition and add to the unit registry
                unit_definition = self._make_cellml_unit(unit_name)
                self.ureg.define(unit_definition)
                self.cellml_defined.add(unit_name)

        # return the defined unit from the registry
        try:
            resolved_unit = self.ureg.parse_expression(unit_name).units
            return resolved_unit
        except pint.UndefinedUnitError:
            raise ValueError('Cannot find the unit with name "%s"' % unit_name)

    def _make_cellml_unit(self, custom_unit_name):
        """Uses the CellML definition for 'unit_name' to construct a Pint unit definition string
        """
        full_unit_expr = []

        # For each of the <unit> elements for this unit definition
        for unit_element in self.cellml_definitions[custom_unit_name]:
            # TODO: handle 'base_units'

            # Source the <unit units="XXX"> from our store
            matched_unit = self.get_quantity(unit_element['units'])

            # Construct a string representing the expression for this <unit>
            expr = str(matched_unit)

            if 'prefix' in unit_element:
                expr = '%s%s' % (unit_element['prefix'], expr)

            if 'exponent' in unit_element:
                expr = '((%s)**%s)' % (expr, unit_element['exponent'])

            # Collect/add this particular <unit> definition
            full_unit_expr.append(expr)

        # Join together all the parts of the unit expression
        full_unit_expr = '*'.join(full_unit_expr)

        # Return Pint definition string
        logger.debug('Unit %s => %s', custom_unit_name, full_unit_expr)
        return '%s = %s' % (custom_unit_name, full_unit_expr)

    @staticmethod
    def one_of_unit(unit):
        """ Returns a quantity of 1 * unit """
        assert isinstance(unit, Unit)
        return 1 * unit

    @staticmethod
    def is_unit_equal(unit1, unit2):
        """ Check whether two Pint Units are equal (converts to quantities if necessary) """
        quantity1 = unit1 if isinstance(unit1, Quantity) else UnitStore.one_of_unit(unit1)
        quantity2 = unit2 if isinstance(unit2, Quantity) else UnitStore.one_of_unit(unit2)
        is_equal = (quantity1.dimensionality == quantity2.dimensionality and
                    math.isclose(quantity1.to(quantity2).magnitude, quantity1.magnitude))
        logger.debug('UnitStore.is_unit_equal(%s, %s) ⟶ %s',
                     quantity1.units, quantity2.units, is_equal)
        return is_equal

    @staticmethod
    def is_quantity_equal(quantity1, quantity2):
        """ Checks whether two instances of Quantity had the same dimensionality and magnitude """
        assert isinstance(quantity1, Quantity)
        assert isinstance(quantity2, Quantity)
        return (quantity1.dimensionality == quantity2.dimensionality
                and quantity1.magnitude == quantity2.magnitude)

    @staticmethod
    def convert_to(unit1, unit2):
        """ Returns a quantity that is the result of converting [one of] unit1 into unit2 """
        assert isinstance(unit1, Unit)
        assert isinstance(unit2, Unit)
        assert unit1.dimensionality == unit2.dimensionality
        quantity1 = unit1 if isinstance(unit1, Quantity) else UnitStore.one_of_unit(unit1)
        return quantity1.to(unit2)

    def summarise_units(self, expr: sympy.Expr):
        """Given a Sympy expression, will get the lambdified string to evaluate units
        """
        to_evaluate = self.printer.doprint(expr)
        to_evaluate = to_evaluate.replace('exp(', 'math.exp(')
        # TODO: get rid of eval
        simplified = eval(to_evaluate, {'u': self.ureg, 'math': math}).units
        logger.debug('UnitStore.summarise_units(%s) ⟶ %s ⟶ %s', expr, to_evaluate, simplified)
        return simplified

    @staticmethod
    def get_conversion_factor(from_unit, to_unit):
        """ Returns the magnitude multiplier required to convert from_unit to to_unit """
        return UnitStore.convert_to(from_unit, to_unit).magnitude


class ExpressionWithUnitPrinter(LambdaPrinter):
    """ A Sympy expression printer that prints the expression with unit information """
    def __init__(self, unit_store: UnitStore = None):
        super().__init__()
        self.unit_store = unit_store

    def __get_dummy_unit(self, expr):
        return self.unit_store.model.dummy_info[expr]['unit']

    def __get_dummy_number(self, expr):
        if 'number' in self.unit_store.model.dummy_info[expr]:
            return self.unit_store.model.dummy_info[expr]['number']
        else:
            return None

    def _print_Dummy(self, expr):
        number = self.__get_dummy_number(expr)
        if number:
            return '%f[%s]' % (number, str(self.__get_dummy_unit(expr)))

        return '%s[%s]' % (expr.name, str(self.__get_dummy_unit(expr)))

    def _print_Derivative(self, expr):
        state_dummy = expr.free_symbols.pop()
        state_unit = self.__get_dummy_unit(state_dummy)
        free_dummy = set(expr.canonical_variables.keys()).pop()
        free_unit = self.__get_dummy_unit(free_dummy)
        return 'Derivative(%s[%s], %s[%s])' % (state_dummy.name,
                                               str(state_unit),
                                               free_dummy.name,
                                               str(free_unit))


class PintUnitPrinter(LambdaPrinter):
    """ A Sympy expression printer that returns Pint unit arithmetic that can be evaluated """
    def __init__(self, unit_store: UnitStore = None):
        super().__init__()
        self.unit_store = unit_store

    def __get_dummy_unit(self, expr):
        return self.unit_store.model.dummy_info[expr]['unit']

    def _print_Dummy(self, expr):
        unit_with_prefix = re.sub(r'\b([a-zA-Z_0-9]+)\b', r'u.\1',
                                  str(self.__get_dummy_unit(expr)))
        return '(1 * (%s))' % unit_with_prefix

    def _print_Derivative(self, expr):
        state_dummy = expr.free_symbols.pop()
        state_unit = self.__get_dummy_unit(state_dummy)
        free_dummy = set(expr.canonical_variables.keys()).pop()
        free_unit = self.__get_dummy_unit(free_dummy)
        return '(1 * u.%s)/(1 * u.%s)' % (state_unit, free_unit)
