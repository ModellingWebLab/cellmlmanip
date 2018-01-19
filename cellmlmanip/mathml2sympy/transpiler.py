"""
Parses Content MathML and returns equivalent SymPy expressions

Content Markup specification: https://www.w3.org/TR/MathML2/chapter4.html
"""
from xml.dom import Node, minidom

import sympy


def parse_string(xml_string):
    """
    Reads MathML content from a string and returns equivalent SymPy expressions
    """
    dom = minidom.parseString(xml_string)
    return parse_dom(dom.childNodes[0])


def parse_dom(math_dom_element):
    """
    Accepts a <math> node of DOM structure and returns equivalent SymPy expressions.
    Note: math_dom_element must point the <math> XmlNode, not the root XmlDocument

    :param math_dom_element: <math> XmlNode object of a MathML DOM structure
    :return: List of SymPy expression(s)
    """
    return transpile(math_dom_element)


def transpile(xml_node):
    """
    Descends the given MathML element node and calls the corresponding handler for child elements.
    Returns the SymPy expression of node
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
                print('Hit text node with text "' + text + '"')
        elif child_node.nodeType == child_node.ELEMENT_NODE:
            # Call the appropriate MathML handler function for this tag
            tag_name = child_node.tagName
            if tag_name in HANDLERS:
                sympy_expressions.append(HANDLERS[tag_name](child_node))
            else:
                # MathML handler function not found for this tag!
                raise NotImplementedError('No handler for element <%s>' % tag_name)
        elif child_node.nodeType not in [Node.COMMENT_NODE, Node.PROCESSING_INSTRUCTION_NODE]:
            raise NotImplementedError('Unknown node type %d' % child_node.nodeType)
    return sympy_expressions


# MATHML ELEMENT HANDLERS ######################################################################

def math_handler(node):
    """
    Descend XML node <math>...</math>
    """
    result = transpile(node)
    return result


# TOKEN ELEMENTS ###############################################################################

def ci_handler(node):
    """
    MathML:  https://www.w3.org/TR/MathML2/chapter4.html#contm.ci
    SymPy: http://docs.sympy.org/latest/modules/core.html#id17
    """
    identifier = node.childNodes[0].data.strip()
    return sympy.Symbol(identifier)


def cn_handler(node):
    """
    MathML: https://www.w3.org/TR/MathML2/chapter4.html#contm.cn
    SymPy: http://docs.sympy.org/latest/modules/core.html#number
    """

    # If this number is using scientific notation
    if 'type' in node.attributes:
        if node.attributes['type'].value == 'e-notation':
            # A real number may also be presented in scientific notation. Such numbers have two
            # parts (a mantissa and an exponent) separated by sep. The first part is a real number,
            # while the second part is an integer exponent indicating a power of the base.
            # For example, 12.3<sep/>5 represents 12.3 times 10^5. The default presentation of
            # this example is 12.3e5.
            if len(node.childNodes) == 3 and node.childNodes[1].tagName == 'sep':
                mantissa = node.childNodes[0].data.strip()
                exponent = int(node.childNodes[2].data.strip())
                return sympy.Float('%se%d' % (mantissa, exponent))
            else:
                raise SyntaxError('Expecting <cn type="e-notation">significand<sep/>exponent</cn>.'
                                  'Got: ' + node.toxml())
        raise NotImplementedError('Unimplemented type attribute for <cn>: '
                                  + node.attributes['type'].value)

    number = float(node.childNodes[0].data.strip())
    return sympy.Number(number)


# BASIC CONTENT ELEMENTS #######################################################################

def apply_handler(node):
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.apply
    """
    result = transpile(node)
    if len(result) > 1:
        expression = result[0](*(result[1:]))
    else:
        expression = result[0]
    return expression


def piecewise_handler(node):
    """
    MathML: https://www.w3.org/TR/MathML2/chapter4.html#contm.piecewise
    SymPy: http://docs.sympy.org/latest/modules/functions/elementary.html#piecewise

    constructor, zero or more <piece>, zero or one <otherwise>
    """
    result = transpile(node)
    return sympy.Piecewise(*result)


def piece_handler(node):
    """
    MathML: https://www.w3.org/TR/MathML2/chapter4.html#contm.piecewise
    Returns a 2-tuple defining an expression and condition
    <piece> element contains exactly two children
    """
    result = transpile(node)
    if len(result) != 2:
        raise ValueError('Need exactly 2 children for <piece>')
    return result[0], result[1]


def otherwise_handler(node):
    """
    MathML: https://www.w3.org/TR/MathML2/chapter4.html#contm.piecewise
    Returns a 2-tuple defining an expression and condition
    """
    result = transpile(node)
    if len(result) != 1:
        raise ValueError('More than 1 child for <otherwise>')
    return result[0], True


# ARITHMETIC, ALGEBRA AND LOGIC ################################################################

def minus_handler(node):
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.minus
    unary arithmetic operator OR binary arithmetic operator

    From http://docs.sympy.org/latest/tutorial/manipulation.html:
    "There is no subtraction class in SymPy. x - y is represented as x + -y, or, more
    completely, x + -1*y, i.e., Add(x, Mul(-1, y))."

    * Negation (-a) is equivalent to sympy.Mul(sympy.S.NegativeOne, a)
    * Subtraction (a - b) is equivalent to sympy.Add(a, sympy.Mul(sympy.S.NegativeOne, b))
    """
    def wrapped_minus(left_operand, right_operand=None):
        if right_operand is None:
            # unary arithmetic operator => negation
            return -left_operand
        # otherwise, binary arithmetic operator => subtraction
        return left_operand - right_operand
    return wrapped_minus


def divide_handler(node):
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.divide
    binary arithmetic operator
    There is no class in SymPy for division. Rather, division is represented by a power of -1.

    Equivalent to sympy.Mul(a, sympy.Pow(b, sympy.S.NegativeOne))
    """
    def wrapped_divide(dividend, divisor):
        return dividend / divisor
    return wrapped_divide


def power_handler(node):
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.power
    binary arithmetic operator
    equivalent to sympy.Pow(a, b)
    """
    def wrapped_power(base, exponent):
        return base ** exponent
    return wrapped_power


def root_handler(node):
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.root
    operator taking qualifiers

    Nasty:
    The root element is used to construct roots. The kind of root to be taken is specified by a
    degree element, which should be given as the second child of the apply element enclosing the
    root element. Thus, square roots correspond to the case where degree contains the value 2, cube
    roots correspond to 3, and so on. If no degree is present, a default value of 2 is used.
    """
    def wrapped_root(first_argument, second_argument=None):
        # if no <degree> given, it's sqrt
        if second_argument is None:
            return sympy.root(first_argument, 2)
        return sympy.root(second_argument, first_argument)
    return wrapped_root


def degree_handler(node):
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.degree
    Meaning of <degree> depends on context! We implement it for order of <bvar> in <diff> and
    the kind of root in <root>
    """
    result = transpile(node)
    if len(result) != 1:
        raise ValueError('Expected single value in <degree> tag.'
                         'Got: ' + node.toxml())
    return result[0]


# CALCULUS AND VECTOR CALCULUS #################################################################

def diff_handler(node):
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.diff
    operator taking qualifiers
    """
    def wrapped_diff(x_symbol, y_symbol, evaluate=False):
        # dx / dy
        y_function = sympy.Function(y_symbol.name)

        # if bound variable element <bvar> contains <degree>, argument x_symbol is a list,
        # otherwise, it is a symbol
        if isinstance(x_symbol, list) and len(x_symbol) == 2:
            bound_variable = x_symbol[0]
            order = int(x_symbol[1])
            return sympy.Derivative(y_function(bound_variable), bound_variable, order,
                                    evaluate=evaluate)

        return sympy.Derivative(y_function(x_symbol), x_symbol, evaluate=evaluate)
    return wrapped_diff


def bvar_handler(node):
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.bvar
    NASTY: bvar element depends on the context it is being used
    In a derivative, it indicates the variable with respect to which a function is being
    differentiated.

    The bound variable <bvar> can also specify degree. In this case, we'll have two elements
    """
    result = transpile(node)
    if len(result) == 1:
        # Bound variable without specifying degree
        return result[0]
    elif len(result) == 2:
        return result
    else:
        raise SyntaxError("Don't know how to handle <bvar> " + node.toxml())


# ELEMENTARY CLASSICAL FUNCTIONS ###############################################################

def log_handler(node):
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.log
    operator taking qualifiers or a unary calculus operator
    """
    def wrapped_log(first_element, second_element=None):
        if second_element is None:
            # if no <logbase> element is present, the base is assumed to be 10
            return sympy.log(first_element, 10)

        # Has <logbase> element, which is the first_element after <log/>
        return sympy.log(second_element, first_element)
    return wrapped_log


def logbase_handler(node):
    """
    Qualifier for <log>

    The log function accepts only the logbase schema. If present, the logbase schema denotes the
    base with respect to which the logarithm is being taken. Otherwise, the log is assumed to be b
    ase 10. When used with log, the logbase schema is expected to contain a single child schema;
    otherwise an error is generated.

    Should be the first element following log, i.e. the second child of the containing apply
    element.
    """
    return transpile(node)[0]


def simple_operator_handler(node):
    """
    This function handles simple MathML <tagName> to sympy.Class operators, where no unique handling
    of tag children etc. is required.
    """
    return getattr(sympy, SIMPLE_MATHML_TO_SYMPY_NAMES[node.tagName])


# END OF MATHML HANDLERS #######################################################################

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

# Mapping MathML tag element names (keys) to appropriate handler for SymPy output (values)
# These tags require explicit handling because they have children or context etc.
HANDLERS = {
    'apply': apply_handler,
    'bvar': bvar_handler,
    'ci': ci_handler,
    'cn': cn_handler,
    'degree': degree_handler,
    'diff': diff_handler,
    'divide': divide_handler,
    'log': log_handler,
    'logbase': logbase_handler,
    'math': math_handler,
    'minus': minus_handler,
    'otherwise': otherwise_handler,
    'piece': piece_handler,
    'piecewise': piecewise_handler,
    'power': power_handler,
    'root': root_handler
}

# Add tags that can be handled by simple_operator_handler
for tagName in SIMPLE_MATHML_TO_SYMPY_NAMES:
    HANDLERS[tagName] = simple_operator_handler
