"""Unit handling for CellML models, using the Pint unit library (replaces previous
Sympy-units implementation
"""
import logging
import math
import numbers
import os
from functools import reduce
from operator import mul
from typing import Dict, List

import pint
import sympy
from pint.converters import ScaleConverter
from pint.definitions import UnitDefinition
from sympy.printing.lambdarepr import LambdaPrinter


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# The full list of supported CellML units
# Taken from https://www.cellml.org/specifications/cellml_1.1/#sec_units
CELLML_UNITS = {
    # Base SI units
    'ampere', 'candela', 'kelvin', 'kilogram', 'meter', 'mole', 'second',

    # Derived SI units
    'becquerel', 'celsius', 'coulomb', 'farad', 'gray', 'henry', 'hertz', 'joule',
    'lumen', 'lux', 'newton', 'ohm', 'pascal', 'radian', 'siemens', 'sievert',
    'steradian', 'tesla', 'volt', 'watt', 'weber',

    'katal',

    # Convenience units
    'dimensionless', 'gram', 'liter',

    # Aliases
    'metre', 'litre',
}

CELLML_UNIT_PREFIXES = {
    'yocto': 1e-24,
    'zepto': 1e-21,
    'atto': 1e-18,
    'femto': 1e-15,
    'pico': 1e-12,
    'nano': 1e-9,
    'micro': 1e-6,
    'milli': 1e-3,
    'centi': 1e-2,
    'deci': 1e-1,
    'deca': 1e+1,
    'deka': 1e+1,
    'hecto': 1e2,
    'kilo': 1e3,
    'mega': 1e6,
    'giga': 1e9,
    'tera': 1e12,
    'peta': 1e15,
    'exa': 1e18,
    'zetta': 1e21,
    'yotta': 1e24
}

TRIG_FUNCTIONS = {
    'arccos', 'arccosh',
    'arccot', 'arccoth',
    'arccsc', 'arccsch',
    'arcsec', 'arcsech',
    'arcsin', 'arcsinh',
    'arctan', 'arctanh',
    'cos', 'cosh',
    'cot', 'coth',
    'csc', 'csch',
    'sec', 'sech',
    'sin', 'sinh',
    'tan', 'tanh',
}


class UnitError(Exception):
    """ base class for unit errors
    """
    pass


class UnitsCannotBeCalculatedError(UnitError):
    """ generic invalid units error
    @param expression -- input expression in which the error occurred
    """

    def __init__(self, expression):
        self.expression = expression
        self.message = 'The units of this expression cannot be calculated.'


class UnexpectedMathUnitsError(UnitError):
    """ error when math is outside the set of MathML expected
    @param expression -- input expression in which the error occurred
    """

    def __init__(self, expression):
        self.expression = expression
        self.message = 'The math used by this expression is not supported.'


class InputArgumentsInvalidUnitsError(UnitError):
    """ error when the arguments to a function have incorrect units
    @param expression -- input expression in which the error occurred
    """

    def __init__(self, expression):
        self.expression = expression
        self.message = 'The arguments to this expression should all have the same units.'


class InputArgumentsMustBeDimensionlessError(UnitError):
    """ error when the arguments to a function have incorrect units when they should be dimensionless
    @param expression -- input expression in which the error occurred
    """

    def __init__(self, expression, position=''):
        self.expression = expression
        if position.__len__() == 0:
            self.message = 'The arguments to this expression should be dimensionless.'
        else:
            self.message = 'The %s argument to this expression should be dimensionless.' % position


class InputArgumentMustBeNumberError(UnitError):
    """ error when the input argument must be a number
    @param expression -- input expression in which the error occurred
    @param position -- the position of the argument with error i.e. first/second
    """

    def __init__(self, expression, position):
        self.expression = expression
        self.message = 'The %s argument to this expression should be a number.' % position


class BooleanUnitsError(UnitError):
    """ error when being asked for units of an expression that will return a boolean
    @param expression -- input expression in which the error occurred
    """

    def __init__(self, expression):
        self.expression = expression
        self.message = 'This expression involves boolean values which do not conform to ' \
                       'unit dimensionality rules.'


class UnitStore(object):
    """Wraps the underlying Pint UnitRegistry to provide unit handling for the model's Sympy
    expressions. Getting and checking units is handled by this class.
    """
    def __init__(self, model):
        """Initialise a UnitStore instance; wraps the unit registry and handles addition of new
        unit definitions"""
        cellml_unit_definition = os.path.join(
            os.path.dirname(__file__), 'cellml_units.txt'
        )
        self.ureg: pint.UnitRegistry = pint.UnitRegistry(cellml_unit_definition)

        # units that are defined and added to the unit registry, on top of default cellml units
        self.custom_defined = set()

        # Keep reference to the underlying model, to look up the 'dummy_info' dictionary
        self.model = model

    def add_custom_unit(self, units_name, unit_attributes):
        """Define a new Pint unit definition to the unit registry
        :param units_name: name of unit to add - cannot be an existing CellML unit
        :param unit_attributes: dict object for unit
            {'multiplier': float, 'units': string, 'exponent': integer, 'prefix': string/integer}
            Not all fields necessary but 'units' must match a unit in Pint registry
        """
        assert units_name not in CELLML_UNITS, 'Cannot redefine CellML unit <%s>' % units_name

        unit_definition = self._make_pint_unit_definition(units_name, unit_attributes)

        if unit_definition is not None:
            assert not self._is_unit_defined(units_name), 'Unit <%s> already exists' % units_name

            # check whether the definition is dimensionless (e.g. a dimensionless scaling factor)
            definition = self.ureg.parse_expression(unit_definition.split('=')[1])
            # if the unit is dimensionless
            if definition.dimensionless:
                # get the scale factor by converting the quantity to dimensionless
                scale_factor = definition.to(self.ureg.dimensionless)
                # dimensionless units can't be created using definition strings
                unit_definition = UnitDefinition(units_name, '', (), ScaleConverter(scale_factor))

            self._define_pint_unit(units_name, unit_definition)

    def add_preferred_custom_unit_name(self, units_name, unit_attributes):
        """Set a prefered name for all equivalent custom units. If it does not exist yet, also define a new Pint unit
        definition to the unit registry.
        PLEASE NOTE: not retrospective, assumes you are not adding new expressions to the model.
        It also does not do any unit conversions, only works on equivalent units with different names.
        :param units_name: the prefered name for this unit and all units in the model that are equivalent
        :param unit_attributes:
        """
        try:
            self.add_custom_unit(units_name, unit_attributes)
        except AssertionError:
            pass  # Unit already exists, but that is not a problem
        for dummy in list(self.model.dummy_metadata.items()):
            if self.is_unit_equal(dummy[1].units, getattr(self.ureg, units_name)):
                dummy[1].units = getattr(self.ureg, units_name)

    def add_base_unit(self, units_name):
        """Define a new base unit in the Pint registry
        :param units_name: string name of unit to add
        """
        assert units_name not in CELLML_UNITS, 'Cannot redefine CellML unit <%s>' % units_name
        assert not self._is_unit_defined(units_name), 'Unit <%s> already exists' % units_name
        self._define_pint_unit(units_name, '{name}=[{name}]'.format_map({'name': units_name}))

    def _is_unit_defined(self, name):
        """ Check whether the unit already exists in Pint registry
        :param name: string name of the unit
        :return: True if exists, else False
        """
        try:
            getattr(self.ureg, name)
            return True
        except pint.UndefinedUnitError:
            return False

    def _define_pint_unit(self, units_name, definition_string_or_instance):
        self.ureg.define(definition_string_or_instance)
        self.custom_defined.add(units_name)

    def get_quantity(self, unit_name: str):
        """Returns a pint.Unit from the UnitRegistry with the given name
        :param unit_name:
        :return: pint.Unit
        """
        try:
            resolved_unit = self.ureg.parse_expression(unit_name).units
            return resolved_unit
        except pint.UndefinedUnitError:
            raise KeyError('Cannot find unit <%s> in unit registry' % unit_name)

    def _make_pint_unit_definition(self, units_name, unit_attributes: List[Dict]):
        """Uses the CellML definition to construct a Pint unit definition string
        """
        full_unit_expr = []

        # For each of the <unit> elements for this unit definition
        for unit_element in unit_attributes:
            # Source the <unit units="XXX"> from our store
            matched_unit = self.get_quantity(unit_element['units'])

            # Construct a string representing the expression for this <unit>
            expr = str(matched_unit)

            # See https://www.cellml.org/specifications/cellml_1.1/#sec_units 5.2.2
            # offset, prefix, exponent, and multiplier

            if 'prefix' in unit_element:
                if unit_element['prefix'] in CELLML_UNIT_PREFIXES:
                    power = CELLML_UNIT_PREFIXES[unit_element['prefix']]
                else:
                    # Assume that prefix is an integer - will CellML validation check?
                    power = '1e%s' % unit_element['prefix']
                expr = '(%s * %s)' % (expr, power)

            if 'exponent' in unit_element:
                expr = '((%s)**%s)' % (expr, unit_element['exponent'])

            if 'multiplier' in unit_element:
                expr = '(%s * %s)' % (unit_element['multiplier'], expr)

            # Collect/add this particular <unit> definition
            full_unit_expr.append(expr)

        # Join together all the parts of the unit expression
        full_unit_expr = '*'.join(full_unit_expr)

        # Return Pint definition string
        logger.debug('Unit %s => %s', units_name, full_unit_expr)
        return '%s=%s' % (units_name, full_unit_expr)

    def is_unit_equal(self, unit1, unit2):
        """Compares two units are equivalent by converting them into base units and comparing the
        resulting units and multiplicative factors"""
        assert isinstance(unit1, self.ureg.Unit)
        assert isinstance(unit2, self.ureg.Unit)
        base1 = self.ureg.get_base_units(unit1)
        base2 = self.ureg.get_base_units(unit2)
        is_equal = math.isclose(base1[0], base2[0]) and base1[1] == base2[1]
        logger.debug('is_unit_equal(%s, %s) ⟶ %s', unit1, unit2, is_equal)
        return is_equal

    def convert_to(self, quantity, unit):
        """Returns a quantity that is the result of converting [one of] unit1 into unit2 """
        assert isinstance(quantity, self.ureg.Quantity)
        assert isinstance(unit, self.ureg.Unit)
        return quantity.to(unit)

    def summarise_units(self, expr: sympy.Expr):
        """Given a Sympy expression, will get the lambdified string to evaluate units
        the call to UnitCalculator:traverse will throw an error if units are bad/cannot be calculated
        """
        unit_calculator = UnitCalculator(self.ureg, self.model.dummy_metadata)
        found = unit_calculator.traverse(expr)

        logger.debug('summarise_units(%s) ⟶ %s', expr, found.units)
        return found.units

    def get_conversion_factor(self, quantity, to_unit):
        """Returns the magnitude multiplier required to convert from_unit to to_unit """
        return self.convert_to(quantity, to_unit).magnitude


class UnitCalculator(object):
    """Evaluates a Sympy expression to determine its units. Note: only supports subset of Sympy
    math"""
    def __init__(self, unit_registry, dummy_metadata):
        """:param unit_registry: instance of Pint UnitRegistry
        :param dummy_metadata: a dictionary providing {dummy: MetaDummy} lookup
        """
        self.ureg: pint.UnitRegistry = unit_registry
        self.dummy_metadata = dummy_metadata

    def _check_unit_of_quantities_equal(self, list_of_quantities):
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

    def traverse(self, expr: sympy.Expr):
        """descends the Sympy expression and performs Pint unit arithmetic on sub-expressions
        :param expr: a Sympy expression
        :returns: the quantity (i.e. expression * unit) of the expression
        :throws: KeyError - if variable not found in metadata
                 UnitsCannotBeCalculatedError - if expression just cannot calculate
                 UnexpectedMathUnitsError - if math is not supported
                 BooleanUnitsError - if math returns booleans
                 InputArgumentsMustBeDimensionlessError - if input arguments should be dimensionless
                 InputArgumentsInvalidUnitsError - if input arguments should have same units
                 InputArgumentMustBeNumberError - if one of input arguments should be a number
        NOTE: Sympy will throw exceptions if the expression is badly formed before we even get here
        """
        # collect the units for each argument of the expression (might be sub-expression)
        quantity_per_arg = []

        # weirdly if you create a symbol for a matrix is does not have the is_Piecewise
        # and so crashes rather than hit the catch at the end that should catch anything we dont support
        if expr.is_Matrix:
            raise UnexpectedMathUnitsError(expr)
        # need to work out units of each child atom of the expression
        # piecewise and derivatives are special cases
        if expr.is_Piecewise:
            # Treat each sub-function an argument (& ignore the conditional statements)
            for arg in expr.args:
                quantity_per_arg.append(self.traverse(arg[0]))
        elif expr.is_Derivative:
            # In a derivative expression the first argument is function to be differentiated
            # the second argument is a tuple of which the first is symbol of the differentiator(?)
            quantity_per_arg.append(self.traverse(expr.args[0]))
            quantity_per_arg.append(self.traverse(expr.args[1][0]))
        else:
            # Otherwise, collect the quantity for each argument of the expression
            for arg in expr.args:
                quantity_per_arg.append(self.traverse(arg))

        # if there was an error determining the quantity of an argument quit now
        # we shouldn't ever get here as an exception should already be thrown
        if None in quantity_per_arg:  # pragma: no cover
            raise UnitsCannotBeCalculatedError(str(expr))

        # Terminal atoms in expressions (Integers and Rationals are used by Sympy itself)
        # I have gone through any flag that sympy might use
        # some of which will be redundant - but if we happen to come across it will throw
        # an exception to tell us a model is very unexpected
        if expr.is_Symbol:
            if expr not in self.dummy_metadata:
                raise KeyError('Metadata entry not found for dummy symbol "%s"' % expr)

            metadata = self.dummy_metadata[expr]

            # is this symbol is a placeholder for a number
            if metadata.is_number:
                out = float(metadata.number) * metadata.units
            else:
                #  if this symbol has an initial value (that is not zero)
                if metadata.initial_value and metadata.initial_value != 0.0:
                    # substitute with the initial value for unit arithmetic
                    out = self.ureg.Quantity(float(metadata.initial_value), metadata.units)
                else:
                    # otherwise, keep the symbol
                    out = self.ureg.Quantity(expr, metadata.units)
            return out
        elif expr == sympy.oo:
            return math.inf * self.ureg.dimensionless

        elif expr == sympy.nan:
            return math.nan * self.ureg.dimensionless

        elif expr.is_Number:
            units = self.ureg.dimensionless
            for dummy, metadata in self.dummy_metadata.items():
                if metadata.number and metadata.number == expr and metadata.units:
                    units = metadata.units
                    break
            if expr.is_Integer:
                out = int(expr) * units
            elif expr.is_Rational:
                # NOTE: can't send back Rational(1,2) * u.dimensionless
                # Used by Sympy for root e.g. sqrt(x) == Pow(x, Rational(1,2))
                out = float(expr) * units
            else:
                out = float(expr) * units
            return out

        elif expr.is_Mul:
            # There are no restrictions on the units of operands
            # The result of this operator has units that are the product of the units on
            # the operands. This product may be simplified according to the rules outlined
            # in Appendix C.3.1.
            out = reduce(mul, quantity_per_arg)
            return out

        # Pow is used by Sympy for exponentiating, roots and division
        elif expr.is_Pow:
            base = quantity_per_arg[0]
            exponent = quantity_per_arg[1]

            # exponent must be dimensionless
            if exponent.units != self.ureg.dimensionless:
                logger.critical('Exponent of Pow is not dimensionless %s', expr)
                raise InputArgumentsMustBeDimensionlessError(str(expr), 'second')

            if not isinstance(exponent.magnitude, (sympy.Number, numbers.Number)):
                logger.critical('Exponent of Pow is not a number (is %s): %s',
                                type(exponent.magnitude).__name__,
                                expr)
                raise InputArgumentMustBeNumberError(str(expr), 'second')

            # if base is dimensionless, return is dimensionless
            if base.units == self.ureg.dimensionless:
                return self.ureg.Quantity(base.magnitude**exponent.magnitude,
                                          self.ureg.dimensionless)

            # base is not dimensionless, raise quantity (magnitude and unit) to power
            return base ** exponent

        elif expr.is_Add:
            # These operators, if applied to more than one operand, require all of their operands
            # to have either equivalent units references, as defined in Appendix C.2.1, or to
            # reference units that have dimensional equivalence, as defined in Appendix C.2.2.

            # If all the Add argument quantities are of the same unit
            if self._check_unit_of_quantities_equal(quantity_per_arg):
                # Return the first quantity as representative unit
                out = quantity_per_arg[0]
                return out

            logger.warning('Add args do not have the same unit: %s', expr)
            raise InputArgumentsInvalidUnitsError(str(expr))

        elif expr.is_Piecewise:
            # If unit of each expression in piecewise is the same
            if self._check_unit_of_quantities_equal(quantity_per_arg):
                # return one of them for the unit arithmetic
                out = quantity_per_arg[0]
                return out

            logger.warning('Piecewise args do not have the same unit.')
            raise InputArgumentsInvalidUnitsError(str(expr))

        elif expr.is_Relational:
            # following discussion with Michael we decided that since a
            # variable in cellml can never have a boolean value then
            # we should not encounter expression that return booleans
            # but cellml spec says these should have same arguments so
            # log this
            if not self._check_unit_of_quantities_equal(quantity_per_arg):
                logger.warning('Relational args do not have the same unit: %s', expr)
            logger.critical('Boolean return: %s', expr)
            raise BooleanUnitsError(str(expr))

        elif expr.is_Boolean:
            # following discussion with Michael we decided that since a
            # variable in cellml can never have a boolean value then
            # we should not encounter expression that return booleans
            raise BooleanUnitsError(str(expr))

        elif expr.is_Function:
            # List of functions that have been checked
            if str(expr.func) not in ['cos', 'acos', 'exp', 'floor', 'log', 'Abs', 'tanh', 'ceiling',
                                      'Pow', 'root']:
                logger.warning('Have not checked unit arithmetic for function %s', expr.func)

            # Handle these functions explicitly
            if expr.func == sympy.Abs:
                return abs(quantity_per_arg[0])
            elif expr.func == sympy.floor:
                # if we're flooring a sympy expression
                if isinstance(quantity_per_arg[0].magnitude, sympy.Expr):
                    return 1 * quantity_per_arg[0].units
                return math.floor(quantity_per_arg[0].magnitude) * quantity_per_arg[0].units
            elif expr.func == sympy.ceiling:
                # if we're ceiling a sympy expression
                if isinstance(quantity_per_arg[0].magnitude, sympy.Expr):
                    return 1 * quantity_per_arg[0].units
                return math.ceil(quantity_per_arg[0].magnitude) * quantity_per_arg[0].units
            elif expr.func == sympy.log:
                if self._is_dimensionless(quantity_per_arg[0]):
                    return 1 * self.ureg.dimensionless
                raise InputArgumentsMustBeDimensionlessError(str(expr))
            elif expr.func == sympy.factorial:
                if self._is_dimensionless(quantity_per_arg[0]):
                    return 1 * self.ureg.dimensionless
                raise InputArgumentsMustBeDimensionlessError(str(expr))
            elif expr.func == sympy.exp:
                # requires operands to have units of dimensionless.
                # result of these has units of dimensionless.
                if self._is_dimensionless(quantity_per_arg[0]):
                    # is the operand is a float
                    if isinstance(quantity_per_arg[0].magnitude, float):
                        # return the exponential of the float as dimensionless
                        return self.ureg.Quantity(math.exp(quantity_per_arg[0].magnitude),
                                                  self.ureg.dimensionless)

                    # magnitude contains an unresolved symbol, we lose it here!
                    # we don't lose it - the unit will be dimensionless - we are just not able to
                    # determine the magnitude of the result
                    return 1 * self.ureg.dimensionless

                logger.critical('Exp operand is not dimensionless: %s', expr)
                raise InputArgumentsMustBeDimensionlessError(str(expr))
            # trig. function on any dimensionless operand is dimensionless
            elif str(expr.func) in TRIG_FUNCTIONS:
                if self._is_dimensionless(quantity_per_arg[0]):
                    return 1 * self.ureg.dimensionless
                raise InputArgumentsMustBeDimensionlessError(str(expr))

            # if the function has exactly one dimensionless argument
            if len(quantity_per_arg) == 1 and self._is_dimensionless(quantity_per_arg[0]):
                # assume the result is dimensionless
                return 1 * self.ureg.dimensionless

        elif expr.is_Derivative:
            out = quantity_per_arg[0] / quantity_per_arg[1]
            return out

        # constants in cellml that are all specified as dimensionless
        elif expr == sympy.pi:
            return math.pi * self.ureg.dimensionless
        elif expr == sympy.E:
            return math.e * self.ureg.dimensionless

        raise UnexpectedMathUnitsError(str(expr))


class ExpressionWithUnitPrinter(LambdaPrinter):
    """Sympy expression printer to print expressions with unit information """
    def __init__(self, symbol_info):
        super().__init__()
        self.symbols = symbol_info

    def _get_dummy_unit(self, expr):
        return self.symbols[expr].units

    def _get_dummy_number(self, expr):
        if self.symbols[expr].is_number:
            return self.symbols[expr].number
        return None

    def _print_Dummy(self, expr):
        number = self._get_dummy_number(expr)
        if number is not None:
            return '%f[%s]' % (number, str(self._get_dummy_unit(expr)))

        return '%s[%s]' % (expr.name, str(self._get_dummy_unit(expr)))

    def _print_Derivative(self, expr):
        state_dummy = expr.free_symbols.pop()
        state_unit = self._get_dummy_unit(state_dummy)
        free_dummy = expr.variables[0]
        free_unit = self._get_dummy_unit(free_dummy)
        return 'Derivative(%s[%s], %s[%s])' % (state_dummy.name,
                                               str(state_unit),
                                               free_dummy.name,
                                               str(free_unit))
