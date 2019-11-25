"""Parses Content MathML and returns equivalent SymPy expressions

Content Markup specification: https://www.w3.org/TR/MathML2/chapter4.html
"""
import logging
from xml.dom import Node, minidom

import sympy


logger = logging.getLogger(__name__)


class Transpiler(object):
    """
    Handles conversion of MathmL to Sympy exprerssions.

    :param symbol_generator: An optional method to create expressions for symbols.
        Must have signature ``f(name) -> sympy.Basic``.
    :param number_generator: An optional method to create expressions for numbers with units.
        Must have signature ``f(value, unit) -> sympy.Basic``.
    """

    def __init__(self, symbol_generator=None, number_generator=None):

        # Create simple lambdas for symbol and number generators
        if symbol_generator is None:
            symbol_generator = lambda x: sympy.Symbol(x)
        if number_generator is None:
            number_generator = lambda x, y: sympy.Float(x)

        # Store symbol and number generators
        self.symbol_generator = symbol_generator
        self.number_generator = number_generator

        # Mapping MathML tag element names (keys) to appropriate handler for SymPy output (values)
        # These tags require explicit handling because they have children or context etc.
        self.handlers = {
            'apply': self._apply_handler,
            'bvar': self._bvar_handler,
            'ci': self._ci_handler,
            'cn': self._cn_handler,
            'degree': self._degree_handler,
            'diff': self._diff_handler,
            'divide': self._divide_handler,
            'log': self._log_handler,
            'logbase': self._logbase_handler,
            'math': self._math_handler,
            'minus': self._minus_handler,
            'otherwise': self._otherwise_handler,
            'piece': self._piece_handler,
            'piecewise': self._piecewise_handler,
            'power': self._power_handler,
            'root': self._root_handler
        }

        # Add tags that can be handled by simple_operator_handler
        for tag_name in SIMPLE_MATHML_TO_SYMPY_NAMES:
            self.handlers[tag_name] = self._simple_operator_handler

    def parse_string(self, xml_string):
        """
        Reads MathML content from a string and returns equivalent SymPy expressions.
        :return: A list of SymPy expressions.
        """
        dom = minidom.parseString(xml_string)
        return self.parse_dom(dom.childNodes[0])

    def parse_dom(self, math_dom_element):
        """Accepts a <math> node of DOM structure and returns equivalent SymPy expressions.

        Note: math_dom_element must point the <math> XmlNode, not the root XmlDocument

        :param math_dom_element: <math> XmlNode object of a MathML DOM structure
        :return: A list of SymPy expressions.
        """
        return self.transpile(math_dom_element)

    def transpile(self, xml_node):
        """Descends the given MathML element node and calls the corresponding handler for child
        elements and returns the SymPy expression of node.
        :param xml_node: a DOM element of parsed MathML
        :return: a list of SymPy expressions
        """
        # Collect the parsed expression(s) (i.e. SymPy output) into list
        sympy_expressions = []

        # For each child element of this DOM node
        for child_node in xml_node.childNodes:
            # If this is a Text node (no child nodes)
            if child_node.nodeType == Node.TEXT_NODE:
                # We should be handling text between tags explicitly in the handler
                # (see cn_handler for an example), show a message
                text = child_node.data.strip()
                if text:
                    logger.warning('Unhandled text node in <%s>: "%s"', child_node.tagName, text)
            elif child_node.nodeType == child_node.ELEMENT_NODE:
                # Call the appropriate MathML handler function for this tag
                tag_name = child_node.tagName
                if tag_name in self.handlers:
                    sympy_expressions.append(self.handlers[tag_name](child_node))
                    logger.debug('Transpiled node %s ⟶ %s',
                                 child_node.toxml(), sympy_expressions[-1])
                else:
                    # MathML handler function not found for this tag!
                    raise NotImplementedError('No handler for element <%s>' % tag_name)
            elif child_node.nodeType not in [Node.COMMENT_NODE, Node.PROCESSING_INSTRUCTION_NODE]:
                raise NotImplementedError('Unknown node type %d' % child_node.nodeType)

        return sympy_expressions

    # MATHML ELEMENT HANDLERS ######################################################################

    def _math_handler(self, node):
        """Descend XML node <math>...</math>
        """
        result = self.transpile(node)
        return result

    # TOKEN ELEMENTS ###############################################################################

    def _ci_handler(self, node):
        """MathML:  https://www.w3.org/TR/MathML2/chapter4.html#contm.ci
        SymPy: http://docs.sympy.org/latest/modules/core.html#id17
        """
        identifier = node.childNodes[0].data.strip()
        return self.symbol_generator(identifier)

    def _cn_handler(self, node):
        """MathML: https://www.w3.org/TR/MathML2/chapter4.html#contm.cn
        SymPy: http://docs.sympy.org/latest/modules/core.html#number
        """

        # If this number is using scientific notation
        if 'type' in node.attributes:
            if node.attributes['type'].value == 'e-notation':
                # A real number may also be presented in scientific notation. Such numbers have two
                # parts (a mantissa and an exponent) separated by sep. The first part is a real
                # number, while the second part is an integer exponent indicating a power of the
                # base.. For example, 12.3<sep/>5 represents 12.3 times 10^5. The default
                # presentation of this example is 12.3e5.
                if len(node.childNodes) == 3 and node.childNodes[1].tagName == 'sep':
                    mantissa = node.childNodes[0].data.strip()
                    exponent = int(node.childNodes[2].data.strip())
                    number = float('%se%d' % (mantissa, exponent))
                else:
                    raise SyntaxError('Expecting '
                                      '<cn type="e-notation">significand<sep/>exponent</cn>.'
                                      'Got: ' + node.toxml())
            else:
                raise NotImplementedError('Unimplemented type attribute for <cn>: '
                                          + node.attributes['type'].value)
        else:
            number = float(node.childNodes[0].data.strip())

        # Get units, if given
        # To-do: We're allowing these to _not_ be set for testing only. Maybe remove this option?
        units = node.attributes.get('cellml:units', None)
        if units is not None:
            units = units.value

        return self.number_generator(number, units)

    # BASIC CONTENT ELEMENTS #######################################################################

    def _apply_handler(self, node):
        """https://www.w3.org/TR/MathML2/chapter4.html#contm.apply
        """
        result = self.transpile(node)

        logger.debug('Result of <apply>:\n\t%s\t⟶\t%s', node.toxml(), result)

        if len(result) > 1:
            expression = result[0](*(result[1:]))
        else:
            expression = result[0]
        return expression

    def _piecewise_handler(self, node):
        """MathML: https://www.w3.org/TR/MathML2/chapter4.html#contm.piecewise
        SymPy: http://docs.sympy.org/latest/modules/functions/elementary.html#piecewise

        constructor, zero or more <piece>, zero or one <otherwise>
        """
        result = self.transpile(node)
        return sympy.Piecewise(*result)

    def _piece_handler(self, node):
        """MathML: https://www.w3.org/TR/MathML2/chapter4.html#contm.piecewise
        Returns a 2-tuple defining an expression and condition
        <piece> element contains exactly two children
        """
        result = self.transpile(node)
        if len(result) != 2:
            raise ValueError('Need exactly 2 children for <piece>')
        return result[0], result[1]

    def _otherwise_handler(self, node):
        """MathML: https://www.w3.org/TR/MathML2/chapter4.html#contm.piecewise
        Returns a 2-tuple defining an expression and condition
        """
        result = self.transpile(node)
        if len(result) != 1:
            raise ValueError('More than 1 child for <otherwise>')
        return result[0], True

    # ARITHMETIC, ALGEBRA AND LOGIC ################################################################

    def _minus_handler(self, node):
        """https://www.w3.org/TR/MathML2/chapter4.html#contm.minus
        unary arithmetic operator OR binary arithmetic operator

        From http://docs.sympy.org/latest/tutorial/manipulation.html:
        "There is no subtraction class in SymPy. x - y is represented as x + -y, or, more
        completely, x + -1*y, i.e., Add(x, Mul(-1, y))."

        * Negation (-a) is equivalent to sympy.Mul(sympy.S.NegativeOne, a)
        * Subtraction (a - b) is equivalent to sympy.Add(a, sympy.Mul(sympy.S.NegativeOne, b))
        """
        def _wrapped_minus(left_operand, right_operand=None):
            if right_operand is None:
                # unary arithmetic operator => negation
                return -left_operand
            # otherwise, binary arithmetic operator => subtraction
            return left_operand - right_operand
        return _wrapped_minus

    def _divide_handler(self, node):
        """https://www.w3.org/TR/MathML2/chapter4.html#contm.divide
        binary arithmetic operator
        There is no class in SymPy for division. Rather, division is represented by a power of -1.

        Equivalent to sympy.Mul(a, sympy.Pow(b, sympy.S.NegativeOne))
        """
        def _wrapped_divide(dividend, divisor):
            return dividend / divisor
        return _wrapped_divide

    def _power_handler(self, node):
        """https://www.w3.org/TR/MathML2/chapter4.html#contm.power
        binary arithmetic operator
        equivalent to sympy.Pow(a, b)
        """
        def _wrapped_power(base, exponent):
            return base ** exponent
        return _wrapped_power

    def _root_handler(self, node):
        """https://www.w3.org/TR/MathML2/chapter4.html#contm.root
        operator taking qualifiers

        Nasty:
        The root element is used to construct roots. The kind of root to be taken is specified by a
        degree element, which should be given as the second child of the apply element enclosing the
        root element. Thus, square roots correspond to the case where degree contains the value 2,
        cube roots correspond to 3, and so on. If no degree is present, a default value of 2 is used
        """
        def _wrapped_root(first_argument, second_argument=None):
            # if no <degree> given, it's sqrt
            if second_argument is None:
                return sympy.root(first_argument, 2)
            return sympy.root(second_argument, first_argument)
        return _wrapped_root

    def _degree_handler(self, node):
        """https://www.w3.org/TR/MathML2/chapter4.html#contm.degree
        Meaning of <degree> depends on context! We implement it for order of <bvar> in <diff> and
        the kind of root in <root>
        """
        result = self.transpile(node)
        if len(result) != 1:
            raise ValueError('Expected single value in <degree> tag.'
                             'Got: ' + node.toxml())
        return result[0]

    # CALCULUS AND VECTOR CALCULUS #################################################################

    def _diff_handler(self, node):
        """https://www.w3.org/TR/MathML2/chapter4.html#contm.diff
        operator taking qualifiers
        """
        def _wrapped_diff(x_symbol, y_symbol, evaluate=False):
            # if bound variable element <bvar> contains <degree>, argument x_symbol is a list,
            # otherwise, it is a symbol
            if isinstance(x_symbol, list) and len(x_symbol) == 2:
                bound_variable = x_symbol[0]
                order = int(x_symbol[1])
                deriv = sympy.Derivative(y_symbol, bound_variable, order, evaluate=evaluate)
            # Otherwise, first degree derivative
            else:
                deriv = sympy.Derivative(y_symbol, x_symbol, evaluate=evaluate)

            return deriv

        return _wrapped_diff

    def _bvar_handler(self, node):
        """https://www.w3.org/TR/MathML2/chapter4.html#contm.bvar
        NASTY: bvar element depends on the context it is being used
        In a derivative, it indicates the variable with respect to which a function is being
        differentiated.

        The bound variable <bvar> can also specify degree. In this case, we'll have two elements
        """
        result = self.transpile(node)
        if len(result) == 1:
            # Bound variable without specifying degree
            return result[0]
        elif len(result) == 2:
            return result
        else:
            raise SyntaxError('Do not know how to handle <bvar> ' + node.toxml())

    # ELEMENTARY CLASSICAL FUNCTIONS ###############################################################

    def _log_handler(self, node):
        """https://www.w3.org/TR/MathML2/chapter4.html#contm.log
        operator taking qualifiers or a unary calculus operator
        """
        def _wrapped_log(first_element, second_element=None):
            if second_element is None:
                # if no <logbase> element is present, the base is assumed to be 10
                return sympy.log(first_element, 10)

            # Has <logbase> element, which is the first_element after <log/>
            return sympy.log(second_element, first_element)
        return _wrapped_log

    def _logbase_handler(self, node):
        """Qualifier for <log>

        The log function accepts only the logbase schema. If present, the logbase schema denotes the
        base with respect to which the logarithm is being taken. Otherwise, the log is assumed to be
        base 10. When used with log, the logbase schema is expected to contain a single child schema
        otherwise an error is generated.

        Should be the first element following log, i.e. the second child of the containing apply
        element.
        """
        return self.transpile(node)[0]

    def _get_nary_relation_callback(self, sympy_relation):
        """Wraps the Sympy binary relation to handle n-ary MathML relations

        :param sympy_relation: handle for binary Sympy relation (Eq, Le, Lt, Ge, Gt)
        :return: callback used by the apply_handler to handle n-ary relations
        """
        def _wrapper_relational(*expressions):
            # If the MathML relation is chaining more than 2 expressions
            if len(expressions) > 2:
                # Convert to multiple Sympy binary relations bundled in an 'And' boolean
                relations = []
                for first, second in zip(expressions[:-1], expressions[1:]):
                    relations.append(sympy_relation(first, second))
                return sympy.And(*relations)
            return sympy_relation(*expressions)
        return _wrapper_relational

    def _simple_operator_handler(self, node):
        """This function handles simple MathML <tagName> to sympy.Class operators, where no unique
        handling of tag children etc. is required.
        """
        tag_name = node.tagName

        handler = getattr(sympy, SIMPLE_MATHML_TO_SYMPY_NAMES[tag_name])

        # Some MathML relations allow chaining but Sympy relations are binary operations
        if tag_name in MATHML_NARY_RELATIONS:
            return self._get_nary_relation_callback(handler)

        return handler


# These MathML tags map directly to Sympy classes and don't require any extra handling
SIMPLE_MATHML_TO_SYMPY_NAMES = {
    'abs': 'Abs',
    'and': 'And',
    'arccos': 'acos',
    'arccosh': 'acosh',
    'arccot': 'acot',
    'arccoth': 'acoth',
    'arccsc': 'acsc',
    'arccsch': 'acsch',
    'arcsec': 'asec',
    'arcsech': 'asech',
    'arcsin': 'asin',
    'arcsinh': 'asinh',
    'arctan': 'atan',
    'arctanh': 'atanh',
    'ceiling': 'ceiling',
    'cos': 'cos',
    'cosh': 'cosh',
    'cot': 'cot',
    'coth': 'coth',
    'csc': 'csc',
    'csch': 'csch',
    'eq': 'Eq',
    'exp': 'exp',
    'exponentiale': 'E',
    'false': 'false',
    'floor': 'floor',
    'geq': 'Ge',
    'gt': 'Gt',
    'infinity': 'oo',
    'leq': 'Le',
    'ln': 'ln',
    'lt': 'Lt',
    'max': 'Max',
    'min': 'Min',
    'neq': 'Ne',
    'not': 'Not',
    'notanumber': 'nan',
    'or': 'Or',
    'pi': 'pi',
    'plus': 'Add',
    'rem': 'Mod',
    'sec': 'sec',
    'sech': 'sech',
    'sin': 'sin',
    'sinh': 'sinh',
    'tan': 'tan',
    'tanh': 'tanh',
    'times': 'Mul',
    'true': 'true',
    'xor': 'Xor',
}

# MathML relation elements that are n-ary operators
MATHML_NARY_RELATIONS = {'eq', 'leq', 'lt', 'geq', 'gt'}
