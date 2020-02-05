"""Unit handling for CellML models, using the Pint unit library (replaces previous Sympy-units implementation).
"""
import logging
import math
import numbers
import os
from functools import reduce
from operator import mul

import pint
import sympy
from pint.converters import ScaleConverter
from pint.definitions import UnitDefinition
from sympy.printing.lambdarepr import LambdaPrinter

from . import model


logger = logging.getLogger(__name__)


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
    """Base class for errors relating to calculating units.
    """
    pass


class UnitsCannotBeCalculatedError(UnitError):
    """Generic invalid units error.

    This will be thrown if the expressions or symbols involved in a calculation cannot be calculated.

    :param expression: input expression in which the error occurred
    """

    def __init__(self, expression):
        self.expression = expression
        self.message = 'The units of this expression cannot be calculated.'


class UnexpectedMathUnitsError(UnitError):
    """Invalid units error thrown when math encountered in an expression is outside the subset of MathML expected.

    :param expression: input expression in which the error occurred
    """

    def __init__(self, expression):
        self.expression = expression
        self.message = 'The math used by this expression is not supported.'


class InputArgumentsInvalidUnitsError(UnitError):
    """Invalid units error thrown when the arguments to a function have incorrect units.

    For example it is incorrect to add 1 [meter] to 2 [seconds].

    :param expression: input expression in which the error occurred
    """

    def __init__(self, expression):
        self.expression = expression
        self.message = 'The arguments to this expression should all have the same units.'


class InputArgumentsMustBeDimensionlessError(UnitError):
    """Invalid units error thrown when the arguments to a function have units when the function expects the arguments to
    be dimensionless.

    For example it is incorrect to use a sine function on an argument with units.

    :param expression: input expression in which the error occurred
    """

    def __init__(self, expression, position=''):
        self.expression = expression
        if position.__len__() == 0:
            self.message = 'The arguments to this expression should be dimensionless.'
        else:
            self.message = 'The %s argument to this expression should be dimensionless.' % position


class InputArgumentMustBeNumberError(UnitError):
    """Invalid unit error thrown when the input argument should be a number.

    For example root(x, y) is invalid unless y is a dimensionless number.

    :param expression: input expression in which the error occurred
    :param position: the position of the argument with error i.e. first/second
    """

    def __init__(self, expression, position):
        self.expression = expression
        self.message = 'The %s argument to this expression should be a number.' % position


class BooleanUnitsError(UnitError):
    """Invalid units error when being asked for units of an expression that will return a boolean.

    For example it is incorrect to use the expression 1 [meter] > 0.5 [seconds].

    :param expression: input expression in which the error occurred
    """

    def __init__(self, expression):
        self.expression = expression
        self.message = 'This expression involves boolean values which do not conform to ' \
                       'unit dimensionality rules.'


class UnitStore(object):
    """
    Maps string names to ``pint.Unit`` and ``pint.Quantity`` objects (providing a unit namespace), and provides unit
    checking and conversion.

    Each ``UnitStore`` object contains the units defined in the CellML specification, and more units can be added by the
    user.

    By default, each ``UnitStore`` maintains its own internal ``pint.UnitRegistry``, so that units from different stores
    cannot be compared or interact with each other. To allow comparison and conversion between unit stores, an existing
    ``UnitStore`` can be passed in at construction time, so that the underlying registry will be shared.

    For example::

        a = UnitStore()
        b = UnitStore()

    creates two separated unit stores, each with their own namespace, and without the option to compare units from ``a``
    to units from ``b``.

    But::

        a = UnitStore()
        b = UnitStore(a)
        c = UnitStore(b)

    creates three unit stores, each with their own namespace, but with a single underlying registry for all three
    stores, allowing e.g. units from ``c`` to be converted to units from ``a``.

    :param store: An existing ``UnitStore`` to share a unit registry with.
    """
    _next_id = 0

    def __init__(self, store=None):

        # Assign this unit store a unique id, to keep its namespace separated from any other unit stores sharing the
        # same registry.
        self._id = 'store' + str(UnitStore._next_id) + '_'
        UnitStore._next_id += 1

        # Get unit registry
        if store is None:
            # Create new unit registry, configured for CellML units.
            self._registry = pint.UnitRegistry(os.path.join(os.path.dirname(__file__), 'data', 'cellml_units.txt'))
        else:
            # Share a registry and calculator with the given unit store
            self._registry = store._registry

        # Names of user defined units.
        self._custom_units = set()

        # Create a unit calculator
        self._calculator = UnitCalculator(self)

        # Expose the Unit class
        # TODO Might not be needed in 0.10 anymore
        self.Unit = self._registry.Unit

    def add_unit(self, name, expression):
        """Adds a unit called ``name`` to the unit store, as defined by the string ``expression``.

        For example::

            add_unit('mm', 'meter / 1000')

        :param name: A string name. Names must be unique and cannot overlap with CellML predefined units.
        :param expression: An expression to define the new unit.
        """
        assert name not in CELLML_UNITS, 'Cannot redefine CellML unit <%s>' % name
        assert not self.is_unit_defined(name), 'Unit <%s> already exists' % name

        # TODO ADD NAME PREFIX

        # Dimensionless units can't be created using a string expression.
        # To test if this is a dimensionless unit, parse the string as a Quantity and check if it's dimensionless
        quantity = self._registry.parse_expression(expression)
        if quantity.units == self._registry.dimensionless:
            definition = UnitDefinition(name, '', (), ScaleConverter(quantity.to(self._registry.dimensionless)))
        else:
            definition = expression

        # Add to registry and list of custom units
        self._registry.define(name + '=' + definition)
        self._custom_units.add(name)

    def add_base_unit(self, name):
        """Add a new base unit.

        :param name: A string name.
        """
        assert name not in CELLML_UNITS, 'Cannot redefine CellML unit <%s>' % name
        assert not self.is_unit_defined(name), 'Unit <%s> already exists' % name

        # TODO ADD NAME PREFIX

        # Add to registry and list of custom units
        self._registry.define(name + '=[' + name + ']')
        self._custom_units.add(name)

    def is_unit_defined(self, name):
        """Check if a unit with the given ``name`` exists.

        :param name: string name of the unit
        :returns: True if exists, else False
        """
        # TODO ADD PREFIX

        try:
            getattr(self._registry, name)
            return True
        except pint.UndefinedUnitError:
            return False

    def get_unit(self, name):
        """Retrieves the ``Unit`` object with the given name.

        :param unit_name: string name of the unit
        :returns: A ``Unit`` object.
        :raises KeyError: If the unit is not defined in this unit store.
        """
        # TODO ADD PREFIX

        try:
            # TODO: Look up in unit mapping, then obtain with getattr(). Don't use parse expression because that accepts
            # anything
            return self._registry.parse_expression(name).units
        except pint.UndefinedUnitError:
            raise KeyError('Cannot find unit <%s> in unit registry.' % name)

    def is_unit_equal(self, unit1, unit2):
        """Compares whether two units are equivalent by converting them into base units and comparing the resulting
        units and multiplicative factors.

        :param unit1: the first Unit object to compare
        :param unit2: the second Unit object to compare
        :returns: ``True`` if units are equal, ``False`` otherwise.
        """
        # TODO CHECK IF THESE ARE REALLY UNITS OR QUANTITIES

        assert isinstance(unit1, self._registry.Unit)
        assert isinstance(unit2, self._registry.Unit)
        base1 = self._registry.get_base_units(unit1)
        base2 = self._registry.get_base_units(unit2)
        is_equal = math.isclose(base1[0], base2[0]) and base1[1] == base2[1]
        logger.debug('is_unit_equal(%s, %s) ⟶ %s', unit1, unit2, is_equal)
        return is_equal

    def convert_to(self, quantity, unit):
        """Returns a quantity that is the result of converting [one of] ``unit1`` into ``unit2``.

        :param quantity: the Unit to be converted, multiplied by '1' to form a Quantity object
        :param unit: Unit object into which the first units should be converted
        :returns: a quantity object with the converted unit and corresponding quantity
        """
        # TODO CHECK IF THESE ARE QUANTITIES OR UNITS

        assert isinstance(quantity, self._registry.Quantity)
        assert isinstance(unit, self._registry.Unit)
        return quantity.to(unit)

    def summarise_units(self, expr):
        """
        Given a Sympy expression, will get the lambdified string to evaluate units.

        Note the internal call to ``UnitCalculator:traverse`` will throw an error if units are bad or cannot be
        calculated.

        :param expr: the Sympy expression on which to evaluate units
        :returns: Unit object representing the units of the epression
        :raises KeyError: if variable not found in metadata
        :raises UnitsCannotBeCalculatedError: if expression just cannot calculate
        :raises UnexpectedMathUnitsError: if math is not supported
        :raises BooleanUnitsError: if math returns booleans
        :raises InputArgumentsMustBeDimensionlessError: if input arguments should be dimensionless
        :raises InputArgumentsInvalidUnitsError: if input arguments should have same units
        :raises InputArgumentMustBeNumberError: if one of input arguments should be a number
        """
        # TODO FIX DOCSTRING. WHAT IS A LAMBDIFIED STIRNG????

        found = self._calculator.traverse(expr)

        logger.debug('summarise_units(%s) ⟶ %s', expr, found.units)
        return found.units

    def get_conversion_factor(self, to_unit, from_unit=None, quantity=None, expression=None):
        """Returns the magnitude multiplier required to convert a unit to the specified unit.

        Note this will work on either a unit, a quantity or an expression, but requires only one of these arguments.

        :param to_unit: Unit object into which the units should be converted
        :param from_unit: the Unit to be converted
        :param quantity: the Unit to be converted, multiplied by '1' to form a Quantity object
        :param expression: an expression from which the Unit is evaluated before conversion
        :returns: the magnitude of the resulting conversion factor
        :raises AssertionError: if no target unit is specified or no source unit is specified
        """
        # TODO CHECK IF THESE ARE UNITS OR STRING OR WHAT

        assert to_unit is not None, 'No unit given as target of conversion; to_unit argument is required'
        assert quantity is not None or from_unit is not None or expression is not None, \
            'No unit given as source of conversion; please use one of from_unit, quantity or expression'
        assert [from_unit, quantity, expression].count(None) == 2, \
            'Multiple target specified; please use only one of from_unit, quantity or expression'

        if from_unit is not None:
            assert isinstance(from_unit, self._registry.Unit), 'from_unit must be of type pint:Unit'
            conversion_factor = self.convert_to(1 * from_unit, to_unit).magnitude
        elif quantity is not None:
            assert isinstance(quantity, self._registry.Quantity), 'quantity must be of type pint:Quantity'
            conversion_factor = self.convert_to(quantity, to_unit).magnitude
        else:
            assert isinstance(expression, sympy.Expr), 'expression must be of type Sympy expression'
            conversion_factor = self.convert_to(1 * self.summarise_units(expression), to_unit).magnitude

        if math.isclose(conversion_factor, 1.0):
            return 1.0
        else:
            return conversion_factor

    def dimensionally_equivalent(self, symbol1, symbol2):
        """Returns whether two expressions, symbol1 and symbol2, are dimensionally_equivalent (same units ignoring a
        conversion factor).

        :param symbol1: the first expression to compare
        :param symbol2: the second expression to compare
        :returns: True if units are equal (regardless of quantity), False otherwise
        """
        # TODO IF THEY ARE EXPRESSIONS THE ARGS SHOULD NOT BE CALLED SYMBOL

        try:
            self.get_conversion_factor(from_unit=self.summarise_units(symbol1),
                                       to_unit=self.summarise_units(symbol2))
            return True
        except pint.errors.DimensionalityError:
            return False


class UnitCalculator(object):
    """
    Evaluates a Sympy expression to determine its units.

    Note: only supports a subset of Sympy math.

    :param unit_store: A :class:`UnitStore`.
    """
    def __init__(self, unit_store):
        self._registry = unit_store._registry

    def _check_unit_of_quantities_equal(self, list_of_quantities):
        """Checks whether all units in a list are equivalent.

        :param list_of_quantities: a list of ``pint.Quantity`` objects.
        :returns: boolean indicating whether all units are equivalent
        """
        if len(list_of_quantities) == 1:
            return True

        def _is_equal(quantity1, quantity2):
            assert isinstance(quantity1, self._registry.Quantity)
            assert isinstance(quantity2, self._registry.Quantity)
            base1 = self._registry.get_base_units(1 * quantity1.units)
            base2 = self._registry.get_base_units(1 * quantity2.units)
            return math.isclose(base1[0], base2[0]) and base1[1] == base2[1]

        list_of_quantities = iter(list_of_quantities)
        first = next(list_of_quantities)
        return all(_is_equal(first, rest) for rest in list_of_quantities)

    def _is_dimensionless(self, quantity):
        return quantity.units.dimensionality == self._registry.dimensionless.dimensionality

    def traverse(self, expr):
        """Descends the Sympy expression and performs Pint unit arithmetic on sub-expressions.

        NOTE: Sympy will raise exceptions if the expression is badly formed.

        :param expr: a Sympy expression
        :returns: the quantity (i.e. magnitude(expression) * unit) of the expression
        :raises KeyError: if variable not found in metadata
        :raises UnitsCannotBeCalculatedError: if expression just cannot calculate
        :raises UnexpectedMathUnitsError: if math is not supported
        :raises BooleanUnitsError: if math returns booleans
        :raises InputArgumentsMustBeDimensionlessError: if input arguments should be dimensionless
        :raises InputArgumentsInvalidUnitsError: if input arguments should have same units
        :raises InputArgumentMustBeNumberError: if one of input arguments should be a number
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

            # is this symbol is a placeholder for a number
            if isinstance(expr, model.NumberDummy):
                return float(expr) * expr.units
            elif isinstance(expr, model.VariableDummy):
                if expr.initial_value and expr.initial_value != 0.0:
                    #  if this symbol has an initial value (that is not zero)
                    # substitute with the initial value for unit arithmetic
                    return self._registry.Quantity(float(expr.initial_value), expr.units)
                else:
                    # otherwise, keep the symbol
                    return self._registry.Quantity(expr, expr.units)
            else:   # pragma: no cover
                # An unexpected type of symbol: print anyway for debugging
                return str(expr)
        elif expr == sympy.oo:
            return math.inf * self._registry.dimensionless

        elif expr == sympy.nan:
            return math.nan * self._registry.dimensionless

        elif expr.is_Number:
            units = self._registry.dimensionless
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
            if exponent.units != self._registry.dimensionless:
                logger.critical('Exponent of Pow is not dimensionless %s', expr)
                raise InputArgumentsMustBeDimensionlessError(str(expr), 'second')

            if not isinstance(exponent.magnitude, (sympy.Number, numbers.Number)):
                logger.critical('Exponent of Pow is not a number (is %s): %s',
                                type(exponent.magnitude).__name__,
                                expr)
                raise InputArgumentMustBeNumberError(str(expr), 'second')

            # if base is dimensionless, return is dimensionless
            if base.units == self._registry.dimensionless:
                return self._registry.Quantity(
                    base.magnitude**exponent.magnitude, self._registry.dimensionless)

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
                    return 1 * self._registry.dimensionless
                raise InputArgumentsMustBeDimensionlessError(str(expr))
            elif expr.func == sympy.factorial:
                if self._is_dimensionless(quantity_per_arg[0]):
                    return 1 * self._registry.dimensionless
                raise InputArgumentsMustBeDimensionlessError(str(expr))
            elif expr.func == sympy.exp:
                # requires operands to have units of dimensionless.
                # result of these has units of dimensionless.
                if self._is_dimensionless(quantity_per_arg[0]):
                    # is the operand is a float
                    if isinstance(quantity_per_arg[0].magnitude, float):
                        # return the exponential of the float as dimensionless
                        return self._registry.Quantity(
                            math.exp(quantity_per_arg[0].magnitude), self._registry.dimensionless)

                    # magnitude contains an unresolved symbol, we lose it here!
                    # we don't lose it - the unit will be dimensionless - we are just not able to
                    # determine the magnitude of the result
                    return 1 * self._registry.dimensionless

                logger.critical('Exp operand is not dimensionless: %s', expr)
                raise InputArgumentsMustBeDimensionlessError(str(expr))
            # trig. function on any dimensionless operand is dimensionless
            elif str(expr.func) in TRIG_FUNCTIONS:
                if self._is_dimensionless(quantity_per_arg[0]):
                    return 1 * self._registry.dimensionless
                raise InputArgumentsMustBeDimensionlessError(str(expr))

            # if the function has exactly one dimensionless argument
            if len(quantity_per_arg) == 1 and self._is_dimensionless(quantity_per_arg[0]):
                # assume the result is dimensionless
                return 1 * self._registry.dimensionless

        elif expr.is_Derivative:
            out = quantity_per_arg[0] / quantity_per_arg[1]
            return out

        # constants in cellml that are all specified as dimensionless
        elif expr == sympy.pi:
            return math.pi * self._registry.dimensionless
        elif expr == sympy.E:
            return math.e * self._registry.dimensionless

        raise UnexpectedMathUnitsError(str(expr))


class ExpressionWithUnitPrinter(LambdaPrinter):
    """Sympy expression printer to print expressions with unit information."""

    def _print_NumberDummy(self, expr):
        return '%f[%s]' % (float(expr), str(expr.units))

    def _print_VariableDummy(self, expr):
        return '%s[%s]' % (expr.name, str(expr.units))

    def _print_Derivative(self, expr):
        state = expr.free_symbols.pop()
        freev = expr.variables[0]
        return 'Derivative(%s[%s], %s[%s])' % (state, state.units, freev, freev.units)

