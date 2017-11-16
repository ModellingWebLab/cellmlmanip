"""
Parses Content MathML and returns equivalent SymPy expressions

Content Markup specification: https://www.w3.org/TR/MathML2/chapter4.html

TODO: Validate MathML input
"""
from xml.dom import minidom, Node

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
            name = child_node.tagName
            if name in HANDLERS:
                # If this tag element itself has children
                if child_node.childNodes:
                    # We want to pass the node to the handler, and it will deal with children
                    sympy_expressions.append(HANDLERS[name](child_node))
                else:
                    # This tag has no children
                    sympy_expressions.append(HANDLERS[name]())
            else:
                # MathML handler function not found for this tag!
                raise NotImplementedError('No handler for element <%s>' % child_node.tagName)
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
    TODO: 'type' attribute?
    """
    identifier = node.childNodes[0].data.strip()
    return sympy.Symbol(identifier)


def cn_handler(node):
    """
    MathML: https://www.w3.org/TR/MathML2/chapter4.html#contm.cn
    SymPy: http://docs.sympy.org/latest/modules/core.html#number
    TODO: 'type' attribute?
    """
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

def plus_handler():
    """
    MathML: https://www.w3.org/TR/MathML2/chapter4.html#contm.plus
    n-ary arithmetic operator
    """
    return sympy.Add


def times_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.times
    n-ary arithmetic operator
    """
    return sympy.Mul


def minus_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.minus
    unary arithmetic operator OR binary arithmetic operator

    From http://docs.sympy.org/latest/tutorial/manipulation.html:
    "There is no subtraction class in SymPy. x - y is represented as x + -y, or, more
    completely, x + -1*y, i.e., Add(x, Mul(-1, y))."

    * Negation (-a) is equivalent to sympy.Mul(sympy.S.NegativeOne, a)
    * Subtraction (a - b) is equivalent to sympy.Add(a, sympy.Mul(sympy.S.NegativeOne, b))
    """
    def __wrapped_minus(left_operand, right_operand=None):
        if right_operand is None:
            # unary arithmetic operator => negation
            return -left_operand
        # otherwise, binary arithmetic operator => subtraction
        return left_operand - right_operand
    return __wrapped_minus


def divide_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.divide
    binary arithmetic operator
    There is no class in SymPy for division. Rather, division is represented by a power of -1.

    Equivalent to sympy.Mul(a, sympy.Pow(b, sympy.S.NegativeOne))
    """
    def __wrapped_divide(dividend, divisor):
        return dividend / divisor
    return __wrapped_divide


def rem_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.rem
    binary arithmetic operator
    """
    return sympy.Mod


def floor_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.floor
    unary operator
    """
    return sympy.floor


def abs_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.abs
    unary arithmetic operator
    """
    return sympy.Abs


def power_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.power
    binary arithmetic operator
    equivalent to sympy.Pow(a, b)
    """
    def __wrapped_power(base, exponent):
        return base ** exponent
    return __wrapped_power


def root_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.root
    operator taking qualifiers
    TODO: implement <degree>
    """
    def __wrapped_root(radicand, degree=None):
        if degree is None:
            # by default, sqrt
            return sympy.root(radicand, 2)
        else:
            raise NotImplementedError
            # return sympy.root(b, a)
    return __wrapped_root


# RELATIONS ####################################################################################

def eq_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.eq
    n-ary operator
    """
    return sympy.Eq


def leq_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.leq
    """
    return sympy.Le


def lt_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.lt
    """
    return sympy.Lt


def geq_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.geq
    n-ary relation
    """
    return sympy.Ge


def gt_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.gt
    n-ary relation
    """
    return sympy.Gt


def and_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.and
    n-ary operator
    """
    return sympy.And


def or_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.or
    n-ary operator
    """
    return sympy.Or


# CALCULUS AND VECTOR CALCULUS #################################################################

def diff_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.diff
    operator taking qualifiers
    """
    def __wrapped_diff(x_symbol, y_symbol, evaluate=False):
        # dx / dy
        y_function = sympy.Function(y_symbol.name)  # given by child element <bvar>
        return sympy.Derivative(y_function(x_symbol), x_symbol, evaluate=evaluate)
    return __wrapped_diff


def bvar_handler(node):
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.bvar
    NASTY: bvar element depends on the context it is being used
    In a derivative, it indicates the variable with respect to which a function is being
    differentiated.
    """
    result = transpile(node)
    if len(result) > 1:
        raise NotImplementedError('multiple <bvar> not implemented')
    return result[0]


# ELEMENTARY CLASSICAL FUNCTIONS ###############################################################

def exp_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.exp
    unary arithmetic operator
    """
    return sympy.exp


def ln_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.ln
    unary calculus operator
    """
    return sympy.ln


def log_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.log
    operator taking qualifiers or a unary calculus operator
    TODO: implement <logbase>
    """
    def __wrapped_log(term, base=None):
        if base is None:
            # if no <logbase> element is present, the base is assumed to be 10
            return sympy.log(term, 10)
        else:
            # return sympy.log(b, a)
            raise NotImplementedError
    return __wrapped_log


def cos_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.trig
    unary trigonometric operator
    """
    return sympy.cos


def tanh_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.trig
    unary trigonometric operator
    """
    return sympy.tanh


def arccos_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.trig
    unary trigonometric operator
    """
    return sympy.acos


# CONSTANT AND SYMBOL ELEMENTS #################################################################

def pi_handler():
    """
    https://www.w3.org/TR/MathML2/chapter4.html#contm.pi
    """
    return sympy.pi


# END OF MATHML HANDLERS #######################################################################

# Mapping MathML tag element names (keys) to appropriate function for SymPy output (values)
HANDLERS = {'abs': abs_handler,
            'and': and_handler,
            'apply': apply_handler,
            'arccos': arccos_handler,
            'bvar': bvar_handler,
            'ci': ci_handler,
            'cn': cn_handler,
            'cos': cos_handler,
            'diff': diff_handler,
            'divide': divide_handler,
            'eq': eq_handler,
            'exp': exp_handler,
            'floor': floor_handler,
            'geq': geq_handler,
            'gt': gt_handler,
            'leq': leq_handler,
            'ln': ln_handler,
            'log': log_handler,
            'lt': lt_handler,
            'math': math_handler,
            'minus': minus_handler,
            'or': or_handler,
            'otherwise': otherwise_handler,
            'pi': pi_handler,
            'piece': piece_handler,
            'piecewise': piecewise_handler,
            'plus': plus_handler,
            'power': power_handler,
            'rem': rem_handler,
            'root': root_handler,
            'tanh': tanh_handler,
            'times': times_handler}
