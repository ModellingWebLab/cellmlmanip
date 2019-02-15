"""
Unit handling for CellML models, using the Pint unit library (replaces previous
Sympy-units implementation
"""
import logging
import math
import re

from functools import reduce
from operator import mul

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
            # TODO: what other attributes can a unit have if it have base_units = 'yes'?
            if 'base_units' in unit_element and unit_element['base_units'] == 'yes':
                return '%s = [%s]' % (custom_unit_name, custom_unit_name)

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

    def post(self, _e):

        def check_unit_list_equal(iterator):
            iterator = iter(iterator)
            try:
                first = next(iterator)
            except StopIteration:
                return True
            return all(first.dimensionality == rest.dimensionality and
                       math.isclose(((1*first.units).to(rest.units)).magnitude, 1.0)
                       for rest in iterator)

        quantity_per_arg = []
        for arg in _e.args:
            quantity_per_arg.append(self.post(arg))

        if _e.is_Dummy:
            if 'number' in self.model.dummy_info[_e]:
                r = (float(self.model.dummy_info[_e]['number']) *
                     self.model.dummy_info[_e]['unit'])
            else:
                info = self.model.find_variable({'sympy.Dummy': _e})
                if 'initial_value' in info[0]:
                    r = self.ureg.Quantity(float(info[0]['initial_value']),
                                           self.model.dummy_info[_e]['unit'])
                else:
                    r = self.ureg.Quantity(_e, self.model.dummy_info[_e]['unit'])
            # print('%s is Dummy -> %s' % (_e, repr(r)))
            return r
        elif _e.is_Pow:
            base = quantity_per_arg[0]
            exponent = quantity_per_arg[1]
            # if both are dimensionless
            if (base.units == self.ureg.dimensionless and
                    exponent.units == self.ureg.dimensionless):
                r = self.ureg.Quantity(base.magnitude**exponent.magnitude, self.ureg.dimensionless)
            elif (base.units != self.ureg.dimensionless and
                  exponent.units == self.ureg.dimensionless):
                r = base ** exponent
            else:
                print('could not Pow', sympy.srepr(_e))
                print('that has', quantity_per_arg)
                print('each type is', str([type(x) for x in _e.args]))
                return None
            # print('%s is Pow(%s) -> %s' % (_e,
            #                                ', '.join([repr(x) for x in quantity_per_arg]),
            #                                repr(r)))
            return r
        elif _e.is_Integer:
            r = int(_e) * self.ureg.dimensionless
            # print('%s is Integer -> %s' % (_e, repr(r)))
            return r
        elif _e.is_Rational:
            # TODO: can't send back Rational(1,2) * u.dimensionless
            r = float(_e) * self.ureg.dimensionless
            # print('%s is Rational -> %s' % (_e, repr(r)))
            return r
        elif _e.is_Mul:
            r = reduce(mul, quantity_per_arg)
            # print('%s is Mul -> %s' % (_e, repr(r)))
            return r
        elif _e.is_Add:
            if check_unit_list_equal(quantity_per_arg):
                r = quantity_per_arg[0]
                # print('%s is Add -> %s' % (_e, repr(r)))
                return r
            else:
                print('All items in Add do not have the same unit', quantity_per_arg)
        elif _e.is_Function:
            if str(_e.func) == 'Abs':
                return abs(quantity_per_arg[0])
            # print('%s is a function %s -> %s' % (_e, _e.func, quantity_per_arg))
            if str(_e.func) not in ['cos', 'acos', 'exp', 'floor']:
                print('!!!Check', _e.func)
            # check that unit os dimensional less
            is_dimensionless = ((quantity_per_arg[0]).to(self.ureg.dimensionless).units
                                == self.ureg.dimensionless)
            if len(quantity_per_arg) == 1 and is_dimensionless:
                return 1 * self.ureg.dimensionless
            else:
                print(type(_e), _e.is_Function, _e.func, _e.args)
                raise RuntimeError('HANDLE %s %s' % (_e, sympy.srepr(_e)))
        elif _e == sympy.pi:
            return math.pi * self.ureg.dimensionless
        elif _e.is_Derivative:
            r = quantity_per_arg[0] / quantity_per_arg[1]
            # print('%s is Derivative -> %s' % (_e, repr(r)))
            return r
        else:
            print(type(_e), _e.is_Function, _e.func, _e.args)
            raise RuntimeError('HANDLE %s %s' % (_e, sympy.srepr(_e)))

    def summarise_units(self, expr: sympy.Expr):
        """Given a Sympy expression, will get the lambdified string to evaluate units
        """
        to_evaluate = self.printer.doprint(expr)
        # TODO: get rid of eval
        try:
            # print('Find units of %s' % expr)
            found = self.post(expr)
            # print('using traversal', found.units)
            simplified = eval(to_evaluate, {'u': self.ureg, 'math': math}).units
            # print(simplified, '=?', found.units)
            conversion = (1*simplified).to(found.units)
            assert conversion.units == found.units and math.isclose(conversion.magnitude, 1.0)
        except Exception as e:
            printer = ExpressionWithUnitPrinter(unit_store=self)
            print('the final unit is', repr(self.post(expr)))
            logger.fatal('Could not summaries units: %s', expr)
            logger.fatal('-> %s', printer.doprint(expr))
            logger.fatal('-> %s', to_evaluate)
            simplified = found.units
        logger.debug('summarise_units(%s) ⟶ %s ⟶ %s', expr, to_evaluate, simplified)
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


_known_functions_math = {
    'acos': 'acos',
    'acosh': 'acosh',
    'asin': 'asin',
    'asinh': 'asinh',
    'atan': 'atan',
    'atan2': 'atan2',
    'atanh': 'atanh',
    'ceiling': 'ceil',
    'cos': 'cos',
    'cosh': 'cosh',
    'erf': 'erf',
    'erfc': 'erfc',
    'exp': 'exp',
    'expm1': 'expm1',
    'factorial': 'factorial',
    'floor': 'floor',
    'gamma': 'gamma',
    'hypot': 'hypot',
    'loggamma': 'lgamma',
    'log': 'log',
    'log10': 'log10',
    'log1p': 'log1p',
    'log2': 'log2',
    'sin': 'sin',
    'sinh': 'sinh',
    'Sqrt': 'sqrt',
    'tan': 'tan',
    'tanh': 'tanh'
}  # Not used from ``math``: [copysign isclose isfinite isinf isnan ldexp frexp pow modf
# radians trunc fmod fsum gcd degrees fabs]
_known_constants_math = {
    'Exp1': 'e',
    'Pi': 'pi',
    # Only in python >= 3.5:
    # 'Infinity': 'inf',
    # 'NaN': 'nan'
}


class PintUnitPrinter(LambdaPrinter):
    """ A Sympy expression printer that returns Pint unit arithmetic that can be evaluated """
    def __init__(self, unit_store: UnitStore = None):
        super().__init__()
        self.unit_store = unit_store
        self.in_pow = False

        for k in _known_functions_math:
            setattr(PintUnitPrinter, '_print_%s' % k, self._print_known_func)

    def _print_known_func(self, expr):
        known = _known_functions_math[expr.__class__.__name__]
        return 'math.{name}({args})'.format(name=known, args=', '.join(map(self._print, expr.args)))

    def _print_Abs(self, expr):
        return 'abs({args})'.format(args=self._print(expr.args[0]))

    def __get_dummy_unit(self, expr):
        # logger.debug('__get_dummy_unit(%s)', expr)
        return self.unit_store.model.dummy_info[expr]['unit']

    def __get_dummy_number(self, expr):
        if 'number' in self.unit_store.model.dummy_info[expr]:
            return self.unit_store.model.dummy_info[expr]['number']
        else:
            return None

    def _print_Dummy(self, expr):
        number = self.__get_dummy_number(expr)
        unit = str(self.__get_dummy_unit(expr))
        unit_with_prefix = re.sub(r'\b([a-zA-Z_0-9]+)\b', r'u.\1', unit)

        if number:
            return '(%f * (%s))' % (number, unit_with_prefix)
        else:
            return '(1*%s)' % unit_with_prefix

    def _print_Derivative(self, expr):
        # logger.debug('_print_Derivative(%s)', expr)
        state_dummy = expr.free_symbols.pop()
        state_unit = self.__get_dummy_unit(state_dummy)
        free_dummy = set(expr.canonical_variables.keys()).pop()
        free_unit = self.__get_dummy_unit(free_dummy)
        return '((1 * u.%s)/(1 * u.%s))' % (state_unit, free_unit)

    def _print_Pow(self, expr, rational=False):
        if expr.exp == 0.5:
            # return '{0}({1})'.format('math.sqrt', self._print(expr.base))
            base = self._print(expr.base)
            return '(({0})**(1/2))'.format(base)
        else:
            base = super(LambdaPrinter, self)._print_Pow(expr)
            return base

    def _print_Pi(self, expr):
        return 'math.pi'
