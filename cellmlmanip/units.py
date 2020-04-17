"""Unit handling for CellML models, using the Pint unit library."""
import logging
import math
import numbers
import os
import re
from functools import reduce
from operator import mul

import pint
import sympy
from pint.converters import ScaleConverter
from pint.definitions import UnitDefinition
from pint.errors import DimensionalityError

from . import model


logger = logging.getLogger(__name__)


# The full list of supported CellML units
# Taken from https://www.cellml.org/specifications/cellml_1.1/#sec_units
_CELLML_UNITS = {
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

_TRIG_FUNCTIONS = {
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


# Regex to find unit names in expressions
# It says (1) not preceded by a number or period, (2) then at least one ascii letter, (3) then numbers or _ allowed too
_WORD = re.compile('(?<![0-9.])[a-zA-Z_]+[a-zA-Z0-9_]*')

# Regex to remove UnitStore prefixes
_STORE_PREFIX = re.compile('(?<![a-zA-Z0-9_])store[0-9]+_')


class UnitError(Exception):
    """Base class for errors relating to calculating units."""
    def __str__(self):
        """Display this error message."""
        return str(self.expression) + ': ' + self.message

    def add_context(self, expression, message):
        """Add extra context to the error message.

        Will append ``('Context: ' + message).format(expression)`` to the error message.

        :param expression: the wider expression within which a unit calculation error occurred
        :param message: further details explaining the context of the error
        """
        self.context_expression = expression
        self.message += ('. Context: ' + message).format(expression)


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
    :param message: optional message with further detail
    """

    def __init__(self, expression, message=''):
        self.expression = expression
        self.message = message or 'The math used by this expression is not supported.'


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


class UnitConversionError(UnitError):
    """Represents failure to convert between incompatible units.

    :param expression: the Sympy expression in which the error occurred
    :param from_unit: the units the expression is in
    :param to_unit: the units we tried to convert to
    """
    def __init__(self, expression, from_units, to_units):
        self.expression = expression
        self.message = 'Cannot convert units from {} to {}'.format(from_units, to_units)


class UnitStore(object):
    """
    Creates and stores units, and has functions for unit conversion.

    Within a UnitStore ``store``, units are represented as ``store.Unit`` objects, while numbers with units are
    represented as ``store.Quantity`` objects. Both classes are inherited from the ``pint`` package.

    Each UnitStore has its own unique namespace for storing units.

    By default, each ``UnitStore`` maintains an internal ``pint.UnitRegistry`` and units from different stores cannot be
    compared or interact with each other. To allow comparison and conversion between unit stores, an existing
    ``UnitStore`` can be passed in at construction time, so that the underlying registry will be shared. (This will
    allow units from different stores to interact, but the unit stores will maintain independent namespaces).

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

    :param store: An existing :class:`UnitStore` to share a unit registry with.
    """
    _next_id = 0

    def __init__(self, store=None):

        # Assign this unit store a unique id, to keep its namespace separated from any other unit stores sharing the
        # same registry.
        self._id = UnitStore._next_id
        UnitStore._next_id += 1
        self._prefix = 'store' + str(self._id) + '_'

        # Get unit registry
        if store is None:
            # Create new unit registry, configured for CellML units.
            self._registry = pint.UnitRegistry(os.path.join(os.path.dirname(__file__), 'data', 'cellml_units.txt'))
        else:
            # Share a registry and calculator with the given unit store
            self._registry = store._registry

        # Names of units known to this unit store (including CellML predefined ones). Stored as the user sees them, so
        # without the prefix to make them unique in the registry.
        self._known_units = set(_CELLML_UNITS)

        # Create a unit calculator
        self._calculator = UnitCalculator(self)

        # Expose the Unit and Quantity classes
        # TODO Might not be needed in 0.10 anymore
        self.Unit = self._registry.Unit
        self.Quantity = self._registry.Quantity

    def add_unit(self, name, expression):
        """Adds a unit called ``name`` to the unit store, as defined by the string ``expression``.

        For example::

            add_unit('mm', 'meter / 1000')

        :param name: A string name. Names must be unique and cannot overlap with CellML predefined units.
        :param expression: An expression to define the new unit.
        :returns: The newly created ``Unit`` object.
        """
        if name in _CELLML_UNITS:
            raise ValueError('Cannot redefine CellML unit <%s>.' % name)
        if name in self._known_units:
            raise ValueError('Cannot redefine unit <%s>.' % name)

        # Add prefix to name
        qname = self._prefix_name(name)

        # Add prefixes inside expression
        expression = _WORD.sub(self._prefix_expression, expression)

        # Dimensionless units can't be created using a string expression.
        # To test if this is a dimensionless unit, parse the string as a Quantity and check if it's dimensionless
        quantity = self._registry.parse_expression(expression)
        if quantity.units == self._registry.dimensionless:
            definition = UnitDefinition(qname, '', (), ScaleConverter(quantity.to(self._registry.dimensionless)))
        else:
            definition = qname + '=' + expression

        # Add to registry
        self._registry.define(definition)

        # Add original name to list of known units
        self._known_units.add(name)

        # Return new unit
        return getattr(self._registry, qname)

    def add_base_unit(self, name):
        """Add a new base unit.

        :param name: A string name.
        """
        if name in _CELLML_UNITS:
            raise ValueError('Cannot redefine CellML unit <%s>.' % name)
        if name in self._known_units:
            raise ValueError('Cannot redefine unit <%s>.' % name)

        # Add prefix to name
        qname = self._prefix_name(name)

        # Add to registry
        self._registry.define(qname + '=[' + qname + ']')

        # Add original name to list of known units
        self._known_units.add(name)

        # Return new unit
        return getattr(self._registry, qname)

    def convert(self, quantity, unit):
        """Converts the given ``quantity`` to be in the given ``unit``.

        :param quantity: a ``Quantity`` to convert.
        :param unit: a ``Unit`` to convert it to.
        :returns: the converted ``Quantity``.
        """
        assert isinstance(quantity, self._registry.Quantity)
        assert isinstance(unit, self._registry.Unit)
        return quantity.to(unit)

    def format(self, unit, base_units=False):
        """
        Returns a string representation of a unit.

        :param unit: The ``Unit`` object to format.
        :param base_units: If this optional argument is set to ``True`` the string will show the base units expansion of
            the given unit.
        """
        if base_units:
            base = self._registry.get_base_units(unit)
            text = str(base[0]) + ' ' + str(base[1])
            # Remove unit store prefixes from any unit store
            return _STORE_PREFIX.sub('', text)
        else:
            return _STORE_PREFIX.sub('', str(unit))

    def get_unit(self, name):
        """Retrieves the ``Unit`` object with the given name.

        :param unit_name: string name of the unit
        :returns: A ``Unit`` object.
        :raises KeyError: If the unit is not defined in this unit store.
        """
        if name not in self._known_units:
            raise KeyError('Unknown unit ' + str(name) + '.')

        return getattr(self._registry, self._prefix_name(name))

    def is_defined(self, name):
        """Check if a unit with the given ``name`` exists."""
        return name in self._known_units

    def is_equivalent(self, unit1, unit2):
        """Tests whether two units are equivalent.

        :param unit1: The first ``Unit`` object to compare.
        :param unit2: The second ``Unit`` object to compare.
        :returns: ``True`` if both units are equivalent, ``False`` otherwise.
        """
        assert isinstance(unit1, self.Unit)
        assert isinstance(unit2, self.Unit)
        base1 = self._registry.get_base_units(unit1)
        base2 = self._registry.get_base_units(unit2)
        is_equal = base1[1] == base2[1] and math.isclose(base1[0], base2[0])
        logger.debug('is_equal(%s, %s) ⟶ %s', unit1, unit2, is_equal)
        return is_equal

    def evaluate_units(self, expr):
        """
        Evaluates and returns the ``Unit`` a Sympy expression is in.

        :param expr: the Sympy expression whose units to evaluate.
        :returns: The calculated ``Unit`` object.
        :raises UnitError: if there are unit errors when evaluating the expression's units.
        """
        found = self._calculator.traverse(expr).units
        logger.debug('evaluate_units(%s) ⟶ %s', expr, found)
        return found

    def get_conversion_factor(self, to_unit, from_unit):
        """Returns the magnitude multiplier required to convert a unit to the specified unit.

        :param to_unit: Unit object into which the units should be converted
        :param from_unit: the Unit to be converted
        :returns: the magnitude of the resulting conversion factor
        """
        assert isinstance(from_unit, self._registry.Unit), 'from_unit must be a unit, not ' + str(from_unit)
        cf = self.convert(1 * from_unit, to_unit).magnitude
        return 1.0 if isinstance(cf, numbers.Number) and math.isclose(cf, 1.0) else cf

    def _prefix_expression(self, match):
        """Takes a regex Match object from _WORD, and adds a prefix (UnitStore id), taking SI prefixes into account."""
        return self._prefix_name(match.group(0))

    def _prefix_name(self, name):
        """Adds a prefix to a unit name."""

        # Note: CellML units don't get prefixes. This is OK because they're the same in any model.
        # It's good because (1) it stops us having to redefine all cellml units, (2) it means 'dimensionless' is still
        # treated in a special way (dimensionless * meter = meter, but special_dimensionless * meter isn't simplified).
        if name in _CELLML_UNITS:
            return name

        return self._prefix + name


class UnitCalculator(object):
    """
    Evaluates a Sympy expression to determine its units.

    Note: only supports a subset of Sympy math.

    :param unit_store: A :class:`UnitStore`.
    """
    def __init__(self, unit_store):
        self._store = unit_store
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

        # Dimensionless is used a lot
        dimensionless = self._store.get_unit('dimensionless')

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
                raise RuntimeError('Unexpected symbol type: ' + str(expr))

        elif expr == sympy.oo:
            return math.inf * dimensionless

        elif expr == sympy.nan:
            return math.nan * dimensionless

        elif expr.is_Number:
            units = dimensionless
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
            if exponent.units != dimensionless:
                logger.critical('Exponent of Pow is not dimensionless %s', expr)
                raise InputArgumentsMustBeDimensionlessError(str(expr), 'second')

            if not isinstance(exponent.magnitude, (sympy.Number, numbers.Number)):
                logger.critical('Exponent of Pow is not a number (is %s): %s',
                                type(exponent.magnitude).__name__,
                                expr)
                raise InputArgumentMustBeNumberError(str(expr), 'second')

            # if base is dimensionless, return is dimensionless
            if base.units == dimensionless:
                return self._registry.Quantity(base.magnitude**exponent.magnitude, dimensionless)

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
                    return 1 * dimensionless
                raise InputArgumentsMustBeDimensionlessError(str(expr))
            elif expr.func == sympy.factorial:
                if self._is_dimensionless(quantity_per_arg[0]):
                    return 1 * dimensionless
                raise InputArgumentsMustBeDimensionlessError(str(expr))
            elif expr.func == sympy.exp:
                # requires operands to have units of dimensionless.
                # result of these has units of dimensionless.
                if self._is_dimensionless(quantity_per_arg[0]):
                    # is the operand is a float
                    if isinstance(quantity_per_arg[0].magnitude, float):
                        # return the exponential of the float as dimensionless
                        return self._registry.Quantity(math.exp(quantity_per_arg[0].magnitude), dimensionless)

                    # magnitude contains an unresolved symbol, we lose it here!
                    # we don't lose it - the unit will be dimensionless - we are just not able to
                    # determine the magnitude of the result
                    return 1 * dimensionless

                logger.critical('Exp operand is not dimensionless: %s', expr)
                raise InputArgumentsMustBeDimensionlessError(str(expr))
            # trig. function on any dimensionless operand is dimensionless
            elif str(expr.func) in _TRIG_FUNCTIONS:
                if self._is_dimensionless(quantity_per_arg[0]):
                    return 1 * dimensionless
                raise InputArgumentsMustBeDimensionlessError(str(expr))

            # if the function has exactly one dimensionless argument
            if len(quantity_per_arg) == 1 and self._is_dimensionless(quantity_per_arg[0]):
                # assume the result is dimensionless
                return 1 * dimensionless

        elif expr.is_Derivative:
            out = quantity_per_arg[0] / quantity_per_arg[1]
            return out

        # constants in cellml that are all specified as dimensionless
        elif expr == sympy.pi:
            return math.pi * dimensionless
        elif expr == sympy.E:
            return math.e * dimensionless

        raise UnexpectedMathUnitsError(str(expr))

    def convert_expression_recursively(self, expr, to_units):
        """Generate a version of the given expression in the requested units.

        Rather than assuming the expression is internally consistent (and hence just wrapping
        in a conversion factor) this will recursively traverse the expression tree and convert
        at each level as needed. Hence if the operands of internal expressions are in dimensionally
        consistent but not equal units, conversions will be applied as needed.

        This method is suitable for use converting the RHS of assignment equations to the units
        desired by the LHS, if the ``Eq`` expression is passed in as ``expr`` and ``to_units`` is
        given as ``None``.

        The conversion strategy for each (sub-)expression depends on the operator:
        - for relational operators, all operands are converted to the units of the first operand
        - for Mul the operands can be in any units, and we convert the result if needed
        - for Add all operands are converted to the desired units (or the units of the first operand
          if no desired units are given)
        - for Pow the exponent must be dimensionless while the operand can be in any units, and we
          convert the result if needed
        - for trig functions, exp, log, etc. the operands have to have dimensionless units
        - for piecewise, the conditions must be dimensionless, and the pieces are set to the desired
          units (or the units of the first piece)
        - for derivatives, numbers, variables, etc. we just convert to the desired units

        :param expr: the Sympy expression to convert
        :param to_units: the desired units of the expression, or ``None`` if we don't care or for
            converting an assignment expression.
        :returns: a Sympy expression in the desired units; the input ``expr`` if no conversion was needed.
        :raises UnitError: if conversion is not possible, using a suitable subclass depending on the exact reason
        """
        try:
            new_expr, was_converted, actual_units = self._convert_expression_recursively(expr, to_units)
        except UnitError as e:
            e.add_context(
                expr,
                'trying to convert "{}" to ' + ('consistent units' if to_units is None else str(to_units)))
            raise
        return new_expr

    def set_lhs_units_from_rhs(self, equation):
        """Set the units of the variable assigned to based on the expression assigned.

        :param equation: a :class:`sympy.Eq` instance with a :class:`VariableDummy` on the left hand side
        :returns: a new version of the equation, but with the RHS using consistent units throughout (with
            conversions if needed) and the ``units`` attribute on the LHS updated to match
        :raises UnitError: if conversion is not possible, using a suitable subclass depending on the exact reason
        """
        lhs = equation.lhs
        assert isinstance(lhs, model.VariableDummy)
        try:
            new_rhs, was_converted, rhs_units = self._convert_expression_recursively(equation.rhs, None)
        except UnitError as e:
            e.add_context(equation, 'trying to set LHS units from those of the RHS in {}')
            raise
        lhs.units = rhs_units
        if was_converted:
            equation = sympy.Eq(lhs, new_rhs)
        return equation

    def _convert_expression_recursively(self, expr, to_units):
        """Helper method for :meth:`convert_expression_recursively` which does the heavy lifting.

        :returns: a tuple ``(new_expr, was_converted, actual_units)``
        """
        print('Trying to convert {} {} to {}'.format(type(expr), expr, to_units))
        was_converted = False  # Tracks whether we needed a units conversion
        dimensionless = self._store.get_unit('dimensionless')

        def maybe_convert_expr(expr, was_converted, to_units, from_units):
            """Helper function that adds a conversion factor if needed.

            :returns: ``(expr or new expr, was_converted, to_units)``
            :raises UnitConversionError: if conversion impossible
            """
            try:
                cf = self._store.get_conversion_factor(to_units, from_units)
            except DimensionalityError:
                raise UnitConversionError(expr, from_units, to_units) from None
            if cf != 1:
                print('Converting {} to {} by factor {}'.format(from_units, to_units, cf))
                was_converted = True
                cf = model.NumberDummy(cf, to_units / from_units)
                expr = cf * expr
            return expr, was_converted, to_units

        def maybe_convert_child(expr, was_converted, to_units):
            """Helper function that does the recursive call and updates ``was_converted``.

            :returns: ``(expr or new expr, was_converted or child was_converted, actual_units)``
            """
            expr, child_was_converted, actual_units = self._convert_expression_recursively(expr, to_units)
            return expr, child_was_converted or was_converted, actual_units

        if expr.is_Matrix:
            # See comment in :meth:`traverse()`!
            raise UnexpectedMathUnitsError(str(expr))
        elif expr.is_Symbol:
            assert isinstance(expr, (model.NumberDummy, model.VariableDummy))
            if to_units is None:
                # Nothing to do, just record the actual units
                actual_units = expr.units
            else:
                expr, was_converted, actual_units = maybe_convert_expr(expr, was_converted, to_units, expr.units)
        elif expr.is_Derivative:
            # Just convert the result if needed
            if (not isinstance(expr.args[0], model.VariableDummy) or
                    len(expr.args) > 2 or
                    expr.args[1][1] > 1 or
                    not isinstance(expr.args[1][0], model.VariableDummy)):
                raise UnexpectedMathUnitsError(
                    expr, 'We only support first order derivatives of single variables wrt time')
            _, was_converted, numerator_units = maybe_convert_child(expr.args[0], was_converted, None)
            _, was_converted, denominator_units = maybe_convert_child(expr.args[1][0], was_converted, None)
            actual_units = numerator_units / denominator_units
            if to_units is not None:
                expr, was_converted, actual_units = maybe_convert_expr(expr, was_converted, to_units, actual_units)
        elif expr.is_Mul:
            # Multiply units of operands, which are allowed to be of any units
            new_args = []
            all_arg_units = []
            for arg in expr.args:
                arg, was_converted, arg_units = maybe_convert_child(arg, was_converted, None)
                new_args.append(arg)
                all_arg_units.append(arg_units)
            actual_units = reduce(mul, all_arg_units)
            if to_units is not None or was_converted:
                # Convert the result if required
                expr = expr.func(*new_args)
                expr, was_converted, actual_units = maybe_convert_expr(expr, was_converted, to_units, actual_units)
        elif expr.is_Pow:
            # Pow is used by Sympy for exponentiating, roots and division
            base, exponent = expr.args
            # Exponent must be dimensionless
            exponent, was_converted, _ = maybe_convert_child(exponent, False, dimensionless)
            try:
                exponent_val = float(exponent)
            except TypeError:
                raise InputArgumentMustBeNumberError(str(expr), 'second') from None
            # Base can be any units, then (try to) convert the result if needed
            base, was_converted, base_units = maybe_convert_child(base, was_converted, None)
            if was_converted:
                expr = expr.func(base, exponent)
            actual_units = base_units ** exponent_val
            if to_units is not None:
                expr, was_converted, actual_units = maybe_convert_expr(expr, was_converted, to_units, actual_units)
        elif expr.is_Add:
            # Convert all arguments to the desired units
            new_args = []
            for arg in expr.args:
                arg, was_converted, actual_units = maybe_convert_child(arg, was_converted, to_units)
                if to_units is None:
                    to_units = actual_units  # Use units of the first argument for subsequent ones
                new_args.append(arg)
            if was_converted and new_args:
                expr = expr.func(*new_args)
        elif expr.is_Relational:
            # The desired units must be dimensionless (or arbitrary)
            if to_units is not None and to_units != dimensionless:
                raise BooleanUnitsError(str(expr) + ' in ' + str(to_units))
            # All arguments converted to the units of the first argument
            arg_units = None
            new_args = []
            for arg in expr.args:
                arg, was_converted, arg_units = maybe_convert_child(arg, was_converted, arg_units)
                new_args.append(arg)
            if was_converted and new_args:
                expr = expr.func(*new_args)
            # Result is dimensionless
            actual_units = dimensionless
        elif expr.is_Piecewise:  # NB: Must come before is_Function!
            # Each condition is dimensionless; each piece gets converted to the desired units
            new_args = []
            for arg in expr.args:
                piece, cond = arg
                new_piece, was_converted, actual_units = maybe_convert_child(piece, was_converted, to_units)
                new_cond, was_converted, _ = maybe_convert_child(cond, was_converted, dimensionless)
                if to_units is None:
                    to_units = actual_units  # Use units of the first piece for subsequent ones
                new_args.append((new_piece, new_cond))
            if was_converted and new_args:
                expr = expr.func(*new_args)
        elif expr.is_Function:
            if expr.func in [sympy.floor, sympy.ceiling, sympy.Abs]:
                # Push any needed conversion to the argument(s)
                actual_units = to_units
            else:
                # All arguments should be dimensionless, as should the whole expr
                if to_units is not None and to_units != dimensionless:
                    raise InputArgumentsMustBeDimensionlessError(str(expr) + ' in ' + str(to_units))
                actual_units = dimensionless
            new_args = []
            for arg in expr.args:
                arg, was_converted, actual_units = maybe_convert_child(arg, was_converted, actual_units)
                new_args.append(arg)
            if was_converted and new_args:
                expr = expr.func(*new_args)
        elif expr.is_number or expr.is_Boolean:  # NB: Must come last
            # Plain numbers can only be dimensionless
            # This will also catch E, pi, oo & nan, true & false
            if to_units is not None and to_units != dimensionless:
                raise InputArgumentsMustBeDimensionlessError(str(expr) + ' in ' + str(to_units))
            actual_units = dimensionless
        else:
            raise UnexpectedMathUnitsError(str(expr))

        print('  Result for {} was converted? {} into {}'.format(expr, was_converted, actual_units))
        return expr, was_converted, actual_units
