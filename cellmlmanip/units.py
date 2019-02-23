"""
Unit handling for CellML models, using the Pint unit library (replaces previous
Sympy-units implementation
"""
import logging
import math
from functools import reduce
from operator import mul

import pint
import sympy
from sympy.printing.lambdarepr import LambdaPrinter


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# The full list of supported CellML units
# Taken from https://www.cellml.org/specifications/cellml_1.1/#sec_units
# Some are not defined by Pint, see comments
CELLML_UNITS = {
    # Base SI units
    'ampere', 'candela', 'kelvin', 'kilogram', 'meter', 'mole', 'second',

    # Derived SI units
    'becquerel', 'celsius', 'coulomb', 'farad', 'gray', 'henry', 'hertz', 'joule',
    'lumen', 'lux', 'newton', 'ohm', 'pascal', 'radian', 'siemens', 'sievert',
    'steradian', 'tesla', 'volt', 'watt', 'weber',

    'katal',  # see __add_custom_units()

    # Convenience units
    'dimensionless', 'gram', 'liter',

    # Aliases
    'metre', 'litre', }

CELLML_UNIT_PREFIXES = {
    'yotta', 'zetta', 'exa', 'peta', 'tera', 'giga', 'mega', 'kilo', 'hecto', 'deka',
    'deci', 'centi', 'milli', 'micro', 'nano', 'pico', 'femto', 'atto', 'zepto', 'yocto'
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
                unit_definition = self._make_pint_unit_definition(unit_name)
                if unit_definition:
                    self.ureg.define(unit_definition)
                self.cellml_defined.add(unit_name)

        # return the defined unit from the registry
        try:
            resolved_unit = self.ureg.parse_expression(unit_name).units
            return resolved_unit
        except pint.UndefinedUnitError:
            raise ValueError('Cannot find the unit with name "%s"' % unit_name)

    def _make_pint_unit_definition(self, cellml_custom_unit):
        """Uses the CellML definition for 'unit_name' to construct a Pint unit definition string
        """
        full_unit_expr = []

        # For each of the <unit> elements for this unit definition
        for unit_element in self.cellml_definitions[cellml_custom_unit]:
            # TODO: what other attributes can a unit have if it have base_units = 'yes'?
            if 'base_units' in unit_element and unit_element['base_units'] == 'yes':
                return '%s = [%s]' % (cellml_custom_unit, cellml_custom_unit)

            # Source the <unit units="XXX"> from our store
            matched_unit = self.get_quantity(unit_element['units'])

            # Construct a string representing the expression for this <unit>
            expr = str(matched_unit)

            # See https://www.cellml.org/specifications/cellml_1.1/#sec_units 5.2.2
            # offset, prefix, exponent, and multiplier

            if 'prefix' in unit_element:
                if unit_element['prefix'] in CELLML_UNIT_PREFIXES:
                    expr = '%s%s' % (unit_element['prefix'], expr)
                else:
                    # Assume that prefix is an integer - will CellML validation check?
                    expr = '(%s * 10**%s)' % (expr, unit_element['prefix'])

            if 'exponent' in unit_element:
                expr = '((%s)**%s)' % (expr, unit_element['exponent'])

            if 'multiplier' in unit_element:
                expr = '(%s * %s)' % (unit_element['multiplier'], expr)

            # Collect/add this particular <unit> definition
            full_unit_expr.append(expr)

        # Join together all the parts of the unit expression
        full_unit_expr = '*'.join(full_unit_expr)

        # to avoid recursion due to pint prefix magic
        if cellml_custom_unit == full_unit_expr:
            return None

        # Return Pint definition string
        logger.debug('Unit %s => %s', cellml_custom_unit, full_unit_expr)
        return '%s = %s' % (cellml_custom_unit, full_unit_expr)

    def one_of_unit(self, unit):
        """ Returns a quantity of 1 * unit """
        assert isinstance(unit, self.ureg.Unit)
        return 1 * unit

    def is_unit_equal(self, unit1, unit2):
        """ Check whether two Pint Units are equal (converts to quantities if necessary) """
        quantity1 = unit1 if isinstance(unit1, self.ureg.Quantity) else self.one_of_unit(unit1)
        quantity2 = unit2 if isinstance(unit2, self.ureg.Quantity) else self.one_of_unit(unit2)
        is_equal = (quantity1.dimensionality == quantity2.dimensionality and
                    math.isclose(quantity1.to(quantity2).magnitude, quantity1.magnitude))
        logger.debug('UnitStore.is_unit_equal(%s, %s) ⟶ %s',
                     quantity1.units, quantity2.units, is_equal)
        return is_equal

    def is_quantity_equal(self, quantity1, quantity2):
        """ Checks whether two instances of Quantity had the same dimensionality and magnitude """
        assert isinstance(quantity1, self.ureg.Quantity)
        assert isinstance(quantity2, self.ureg.Quantity)
        return (quantity1.dimensionality == quantity2.dimensionality
                and quantity1.magnitude == quantity2.magnitude)

    def convert_to(self, unit1, unit2):
        """ Returns a quantity that is the result of converting [one of] unit1 into unit2 """
        assert isinstance(unit1, self.ureg.Unit)
        assert isinstance(unit2, self.ureg.Unit)
        assert unit1.dimensionality == unit2.dimensionality
        quantity1 = unit1 if isinstance(unit1, self.ureg.Quantity) else self.one_of_unit(unit1)
        return quantity1.to(unit2)

    def summarise_units(self, expr: sympy.Expr):
        """Given a Sympy expression, will get the lambdified string to evaluate units
        """

        # TODO: move this!
        # get all the symbols in this expression
        symbols = expr.atoms(sympy.Symbol)

        # collect the initial values, if specified, for each symbol
        subs = {}
        for symbol in symbols:
            symbol_info = self.model.find_variable({'sympy.Dummy': symbol})
            if len(symbol_info) == 1:
                if symbol_info[0].get('initial_value', None):
                    subs[symbol] = float(symbol_info[0]['initial_value'])
            else:
                # dummy placeholders for numbers don't have variable information
                assert 'number' in self.model.dummy_info[symbol]

        unit_calculator = UnitCalculator(self.ureg, self.model.dummy_info, subs)

        try:
            # print('Find units of %s' % expr)
            found = unit_calculator.traverse(expr)
            logger.debug('summarise_units(%s) ⟶ %s', expr, found.units)
            # print('using traversal', found.units)
            # simplified = eval(to_evaluate, {'u': self.ureg, 'math': math}).units
            # print(simplified, '=?', found.units)
            # conversion = (1*simplified).to(found.units)
            # assert conversion.units == found.units and math.isclose(conversion.magnitude, 1.0)
        except Exception as e:
            printer = ExpressionWithUnitPrinter(symbol_info=self.model.dummy_info)
            print('The final unit is', repr(unit_calculator.traverse(expr)))
            logger.fatal('Could not summaries units: %s', expr)
            logger.fatal('-> %s', printer.doprint(expr))
            return None
        return found.units

    def get_conversion_factor(self, from_unit, to_unit):
        """ Returns the magnitude multiplier required to convert from_unit to to_unit """
        return self.convert_to(from_unit, to_unit).magnitude


class UnitCalculator(object):
    """Evaluates a Sympy expression to determine its units. Note: only supports subset of Sympy
    math"""
    def __init__(self, unit_registry, symbol_info, symbol_subs):
        """
        :param unit_registry: instance of Pint UnitRegistry
        :param symbol_info: dictionary mapping symbol to 'unit' and, if applicable, 'number'
        :param symbol_subs: dictionary mapping symbol to its numerical replacement
        """
        self.ureg = unit_registry
        self.symbols = symbol_info
        self.subs = symbol_subs

    def _check_unit_list_equal(self, list_of_quantities):
        """checks whether all units in list are equivalent
        :param list_of_quantities: a list of Pint quantities
        :return: boolean indicating whether all units are equivalent
        """
        if len(list_of_quantities) == 1:
            return True

        list_of_quantities = iter(list_of_quantities)
        first = next(list_of_quantities)

        def _is_equal(quantity1, quantity2):
            assert isinstance(quantity1, self.ureg.Quantity)
            assert isinstance(quantity2, self.ureg.Quantity)
            base1 = self.ureg.get_base_units(1 * quantity1.units)
            base2 = self.ureg.get_base_units(1 * quantity2.units)
            return math.isclose(base1[0], base2[0]) and base1[1] == base2[1]

        return all(_is_equal(first, rest) for rest in list_of_quantities)

    def _is_dimensionless(self, quantity):
        return quantity.units.dimensionality == self.ureg.dimensionless.dimensionality

    def traverse(self, expr):
        """descends the Sympy expression and performs Pint unit arithmetic on sub-expressions
        :param expr: a Sympy expression
        """

        # collect the units for each sub-expression in the expression
        quantity_per_arg = []
        for arg in expr.args:
            quantity_per_arg.append(self.traverse(arg))

        if expr.is_Symbol:
            # is this symbol is a placeholder for a number
            if 'number' in self.symbols[expr]:
                r = float(self.symbols[expr]['number']) * self.symbols[expr]['unit']
            else:
                # if we need to replace this symbol for evaluating units
                if expr in self.subs:
                    r = self.ureg.Quantity(self.subs[expr], self.symbols[expr]['unit'])
                else:
                    # otherwise, straightforward Quantity
                    r = self.ureg.Quantity(expr, self.symbols[expr]['unit'])
            return r
        elif expr.is_Pow:
            base = quantity_per_arg[0]
            exponent = quantity_per_arg[1]
            # if both base and exponent are dimensionless units
            if base.units == exponent.units == self.ureg.dimensionless:
                r = self.ureg.Quantity(base.magnitude**exponent.magnitude, self.ureg.dimensionless)
            # if the base unit is not dimensionless, we exponentiate the unit
            elif (base.units != self.ureg.dimensionless and
                  exponent.units == self.ureg.dimensionless):
                r = base ** exponent
            else:
                raise NotImplemented
            return r
        elif expr.is_Integer:
            r = int(expr) * self.ureg.dimensionless
            return r
        elif expr.is_Rational:
            # TODO: can't send back Rational(1,2) * u.dimensionless
            r = float(expr) * self.ureg.dimensionless
            return r
        elif expr.is_Mul:
            r = reduce(mul, quantity_per_arg)
            return r
        elif expr.is_Add:
            if self._check_unit_list_equal(quantity_per_arg):
                r = quantity_per_arg[0]
                return r
            else:
                logger.warning('Add args do not have the same unit. In base units:')
                for x in quantity_per_arg:
                    logger.warning('%s -> %s' % (x.units, self.ureg.get_base_units(x.units)))
                return None
        elif expr.is_Function:
            # List of functions I've checked
            if str(expr.func) not in ['cos', 'acos', 'exp', 'floor', 'log', 'Abs', 'tanh']:
                logger.warning('Have not check unit arithmetic for function %s', expr.func)
            if expr.func == sympy.Abs:
                return abs(quantity_per_arg[0])
            elif expr.func == sympy.log:
                if self._is_dimensionless(quantity_per_arg[0]):
                    return 1 * self.ureg.dimensionless
                else:
                    raise RuntimeError('log args not dimensionless (%s)',
                                       [x.units for x in quantity_per_arg])

            # if this function has exactly one dimensionless argument
            if len(quantity_per_arg) == 1 and self._is_dimensionless(quantity_per_arg[0]):
                # assume the result is dimensionless!
                return 1 * self.ureg.dimensionless
            else:
                print(expr.func, expr.args, quantity_per_arg)
                raise RuntimeError(' %s %s' % (expr, sympy.srepr(expr)))
        elif expr == sympy.pi:
            return math.pi * self.ureg.dimensionless
        elif expr.is_Derivative:
            r = quantity_per_arg[0] / quantity_per_arg[1]
            return r
        # TODO: Handle _e.is_Piecewise here and remove from model.check_left_right_units_equal
        else:
            print(expr.func, expr.args, quantity_per_arg)
            raise RuntimeError('HANDLE %s %s' % (expr, sympy.srepr(expr)))


class ExpressionWithUnitPrinter(LambdaPrinter):
    """ A Sympy expression printer that prints the expression with unit information """
    def __init__(self, symbol_info=None):
        super().__init__()
        if symbol_info is None:
            symbol_info = dict()
        self.symbols = symbol_info

    def __get_dummy_unit(self, expr):
        return self.symbols[expr]['unit']

    def __get_dummy_number(self, expr):
        if 'number' in self.symbols[expr]:
            return self.symbols[expr]['number']
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
