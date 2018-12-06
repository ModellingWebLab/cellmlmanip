import logging
from typing import List

import pint
import pint.unit
import sympy
from sympy.physics import units as units

logging.basicConfig(level=logging.DEBUG)

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


class QuantityStorePints(object):

    def __init__(self, cellml_def=None):
        """Initialise a QuantityStore instance
        :param cellml_def: a dictionary of <units> definitions from the CellML model. See parser
        for format, essentially: {'name of unit': { [ unit attributes ], [ unit attributes ] } }
        """
        # Initialise the unit registry
        # TODO: create Pint unit definition file for CellML
        self.ureg = pint.UnitRegistry()

        # Add default CellML units not provided by Pint
        self.__add_undefined_units()

        # Hold on to custom unit definitions
        self.cellml_definitions = cellml_def if cellml_def else {}
        self.cellml_defined = set()

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
        logging.debug('Unit %s => %s', custom_unit_name, full_unit_expr)
        return '%s = %s' % (custom_unit_name, full_unit_expr)

    @staticmethod
    def get_conversion_factor(from_unit, to_unit):
        assert isinstance(from_unit, pint.unit._Unit)
        assert isinstance(to_unit, pint.unit._Unit)

        assert from_unit.dimensionality == to_unit.dimensionality

        from_quantity = 1 * from_unit
        to_quantity = 1 * to_unit

        return from_quantity.to(to_quantity).magnitude


class QuantityStore(object):
    """Holds sympy.physics.unit.Quantity objects for the model. Can be initialised with <units>
    definitions from the CellML <model>. The get_quantity() methods will find the right
    quantity to return, given the name of the unit (from the internal store, Sympy built-in, or
    constructed anew)
    """

    # Aliases for units to their Sympy equivalent
    UNIT_ALIASES = {
        'litre': 'liter',
        'metre': 'meter',
    }

    # The full list of supported CellML units
    # Taken from https://www.cellml.org/specifications/cellml_1.1/#sec_units
    # Some are not defined by Sympy, see comments
    CELLML_UNITS = [
        # Base SI units
        'ampere',
        'candela',
        'kelvin',
        'kilogram',
        'meter',
        'metre',  # in self.UNIT_ALIASES
        'mole',
        'second',

        # Derived SI units
        'becquerel',  # see __add_custom_units()
        # 'celsius',  # Not supported in CellML 2.0 nor Sympy
        'coulomb',
        'farad',
        'gray',  # see __add_custom_units()
        'henry',
        'hertz',
        'joule',
        'katal',  # see __add_custom_units()
        'lumen',  # see __add_custom_units()
        'lux',
        'newton',
        'ohm',
        'pascal',
        'radian',
        'siemens',
        'sievert',  # see __add_custom_units()
        'steradian',
        'tesla',
        'volt',
        'watt',
        'weber',

        # Convenience units
        'dimensionless',  # see __add_custom_units()
        'gram',
        'liter',
        'litre',  # in self.UNIT_ALIASES
    ]

    def __init__(self, cellml_def=None):
        """Initialise a QuantityStore instance
        :param cellml_def: a dictionary of <units> definitions from the CellML model. See parser
        for format, essentially: {'name of unit': { [ unit attributes ], [ unit attributes ] } }
        """
        # Initialise the store
        self.store = {}

        # Add required quantities not provided by Sympy
        self.__add_custom_units()

        self.cellml_definitions = cellml_def if cellml_def else {}
        self.sympify_context = {}

        # Setup the context with Sympy unit definitions for sympify
        from sympy.core.compatibility import exec_
        exec_('from sympy.physics.units import *', self.sympify_context)

    def __add_custom_units(self):
        """Adds custom Sympy dimensions and quantities that aren't provided by default but
        required by the CellML specification
        """
        self.store['dimensionless'] = units.Quantity('dimensionless',
                                                     units.Dimension(1),
                                                     sympy.Rational(1, 1))

        self.store['becquerel'] = units.Quantity('becquerel', units.Dimension('activity'), 1, 'Bq')

        # Taken from https://github.com/sympy/sympy/pull/13658
        # gray is a J/kg physical quantity
        self.store['gray'] = units.Quantity('gray',
                                            units.energy/units.mass,
                                            units.meter**2/units.second**2)

        self.store['katal'] = units.Quantity('katal',
                                             units.amount_of_substance/units.time,
                                             units.mol/units.second)

        self.store['lumen'] = units.Quantity('lumen',
                                             units.luminous_intensity * units.steradian.dimension,
                                             units.candela * units.steradian)

        # See https://en.wikipedia.org/wiki/Sievert#Definition for relationship with gray
        # sievert is J/kg biological effect
        self.store['sievert'] = units.Quantity('sievert',
                                               units.energy/units.mass,
                                               units.meter**2/units.second**2)

    def get_quantity(self, unit_name):
        """Given the name of the unit, this will either (i) return a Quantity from the internal
        store as its already been resolved previously (ii) return the Quantity from Sympy if it a
        built-in name or (iii) construct and return a new Quantity object using the
        <units><unit></units> definition in the CellML <model>
        """
        # Correct any aliases (e.g. british -> american spelling)
        unit_name = QuantityStore.UNIT_ALIASES.get(unit_name, unit_name)

        # If we've already sourced the quantity for this name
        if unit_name in self.store:
            return self.store[unit_name]

        # If this unit name is defined in the CellML model
        if unit_name in self.cellml_definitions:
            # Make a Sympy Quantity object for this unit definition
            quantity = self._make_cellml_quantity(unit_name)
            self.store[unit_name] = quantity
            return self.store[unit_name]

        # If this unit name is part of the CellML spec and available as built-in Sympy quantity
        if unit_name in QuantityStore.CELLML_UNITS and hasattr(units, unit_name):
            self.store[unit_name] = getattr(units, unit_name)
            return self.store[unit_name]

        raise RuntimeError('Cannot find the unit with name "%s"' % unit_name)

    def summarise_units(self, expr: sympy.Expr):
        """Given a Sympy expression, will get all the Quantity objects in the expression and
        collect them together to give a single Sympy expression of the units
        """
        # remove all free_symbols TODO: why doesn't this work??
        # replace_symbols = {k: 1 for k in expr.free_symbols}
        # with sympy.evaluate(False):
        #    expr = expr.subs(replace_symbols, simultaneous=True)

        # If this expression is a Quantity itself
        if isinstance(expr, units.Quantity):
            # Return it as it is
            return expr

        # Don't descend into Derivative expression (there are no units within!)
        if expr.is_Derivative:
            return None

        if expr.is_Number:
            return expr

        # If expr is an exponential, the result is always dimensionless (don't descend)
        if isinstance(expr, sympy.exp):
            return self.get_quantity('dimensionless')

        # If expression is of the form dimensionless**n then -> dimensionless
        if (isinstance(expr, sympy.Pow)
                and isinstance(expr.args[0], units.Quantity)
                and expr.args[0] == self.get_quantity('dimensionless')):
            return self.get_quantity('dimensionless')

        # remove negative coefficients (prevents cases where `unit - unit = 0`)
        def __no_neg_coeffs(expression):
            if expression.is_Symbol:
                return expression
            elif expression.is_Number:
                return expression
            else:
                coefficients = {k: v for k, v in expression.as_coefficients_dict().items() if v < 0}
                for c in coefficients.keys():
                    expression = expression.subs(c, c * -1)
                return expression.func(*[__no_neg_coeffs(x) for x in expression.args])
        expr = __no_neg_coeffs(expr)

        # Units are always part of a Multiplicative expression
        # TODO: check if all other units are always multiplicative!
        if expr.is_Mul:
            # We only keep quantities
            keep: List[units.Quantity] = []
            # For each of the multiplication arguments that contain Quantity atoms
            for arg in [x for x in expr.args if x.atoms(units.Quantity)]:
                # Descend into the argument and get the Quantities within
                should_keep = self.summarise_units(arg)
                # Keep anything we find (check - we always should!)
                if should_keep:
                    keep.append(should_keep)
            # Perform the same mathematical operation, but on units only (no symbols or numbers)
            return expr.func(*keep)

        # if there are not quantities in this expression, we can't continue
        # (or we should handle this scenario above)
        if not expr.atoms(units.Quantity):
            raise ValueError('No quantities to summarise in expression %s' % expr)

        # Otherwise, descend into the expression tree
        return expr.func(*[self.summarise_units(x) for x in expr.args])

    def is_equal(self, quantity_1, quantity_2):
        """Converts one quantity into another to see if they are equal

        :param quantity_1: a Sympy Quantity instance
        :param quantity_2: a Sympy Quantity instance
        """
        conversion = units.convert_to(quantity_1, quantity_2)
        return conversion == quantity_2

    @staticmethod
    def get_conversion_factor(from_unit, to_unit):
        """Returns the multiplier to convert one unit to another (of the same dimension)"""
        return units.convert_to(from_unit, to_unit).n()

    def _sympify(self, string):
        """Takes a string containing a Sympy expression and evaluates it using the required context
        for handling our units
        """
        # logging.info('sympy.sympify(%s)', string)
        return sympy.sympify(string, locals=self.sympify_context)

    def _make_cellml_quantity(self, name):
        """Will use the CellML unit definition and construct a new Quantity object for that unit
        """
        full_unit_expr = []
        full_dimension = []

        # For each of the <unit> elements for this unit definition
        for unit_element in self.cellml_definitions[name]:
            # if this unit is a new base unit defined by the cellml model
            if 'base_units' in unit_element and unit_element['base_units'] == 'yes':
                return units.Quantity(name, units.Dimension(1), sympy.Rational(1, 1))

            # Source the <unit units="XXX"> from our store
            unit_as_quantity = self.get_quantity(unit_element['units'])

            # Add this unit to the sympify context if necessary
            if unit_element['units'] not in self.sympify_context:
                self.sympify_context[unit_element['units']] = unit_as_quantity

            # Construct a string representing the expression and dimensions for this <unit>
            expr = str(unit_as_quantity.name)
            dimension = str(unit_as_quantity.args[1].args[0])

            if 'prefix' in unit_element:
                expr = '(%s * %s)' % (unit_element['prefix'], expr)

            if 'exponent' in unit_element:
                expr = '((%s)**%s)' % (expr, unit_element['exponent'])
                # Exponent of the <unit> also affects the dimension
                dimension = '(%s)**%s' % (dimension, unit_element['exponent'])

            # Collect/add this particular <unit> definition
            full_unit_expr.append(expr)
            full_dimension.append('(%s)' % dimension)

        # Have sympy evaluate the collected unit dimension and expression
        full_dimension = self._sympify('*'.join(full_dimension))
        full_unit_expr = self._sympify('*'.join(full_unit_expr))

        # Create and return the Quantity object for this CellML <units>
        quantity = units.Quantity(name, full_dimension, full_unit_expr)
        logging.debug('%s=Quantity("%s", %s, %s)',
                      name, name, full_dimension.args[0], full_unit_expr)
        return quantity

    def simplify_units_until_no_change(self, expr):
        """Simplifies the units of an expression until they cannot be simplified any further"""
        unsimplified_expr = expr
        while True:
            simplified_expr = self.summarise_units(unsimplified_expr)
            if unsimplified_expr == simplified_expr:
                break
            unsimplified_expr = simplified_expr
        return simplified_expr