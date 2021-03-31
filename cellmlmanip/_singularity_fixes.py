"""
``remove_fixable_singularities`` specifies a method to remove fixable singularities from a model
"""
from math import isclose

import networkx as nx
from sympy.sets.conditionset import ConditionSet
from sympy.sets.sets import Complement
from sympy import (
    Abs,
    Add,
    Eq,
    Float,
    I,
    Mul,
    Piecewise,
    Pow,
    Wild,
    exp,
    log,
    solveset,
    solve,
    I
)

from .model import FLOAT_PRECISION, Quantity

ONE = Quantity(1.0, 'dimensionless')


def _generate_piecewise(expr, V, sp, Vmin, Vmax):
    """
    Generates a new (piecewise) expression based on expr with a linear interpolation when V is between Vmin and Vmax

    The returned expression between Vmin and Vmax is::
    f(Vmin) + (V - Vmin) / (Vmax - Vmin) * (f(Vmax) - f(Vmin))
where ``f`` is ``expr``

    :param expr: The expression to turn into a piecewise
    :param V: The voltage variable.
    :param sp: The value of the singularity point: the value of V for which U == 0.
    :param Vmin: The value of the lower bound: the value of V for which U - U_offset == 0.
    :param Vmax: The value of the upper bound: the value of V for which U + U_offset == 0.
    :return: expr with a singularity exception around sp, if sp is not None else return the original expr.
    """
    f_Vmin = expr.xreplace({V: Vmin})
    f_Vmax = expr.xreplace({V: Vmax})
    return Piecewise((f_Vmin + ((V - Vmin) / (Vmax - Vmin)) * (f_Vmax - f_Vmin),
                      Abs(V - sp) < Abs((Vmax - Vmin) / 2)), (expr, True))


def _float_dummies(expr):
    """Turns floats back into Quantity dummies to be in line with the rest of the sympy equations"""
    if expr is None:
        return None
    return expr.xreplace({f: Quantity(f, 'dimensionless') for f in expr.atoms(Float)})


def _is_negative_power(expr):
    try:
        return isinstance(expr, Pow) and bool(expr.args[1] < 0)
    except TypeError:  # if this is a power with variables still in it we can't determine if it's negative
        pass
    return False


def _get_singularity(expr, V, U_offset, exp_function):
    """
    Finds singularities in equations of the form:
       - ``U / (exp(U) - 1.0)``
       - ``U / (1.0 - exp(U))``
       - ``(exp(U) - 1.0) / U``
       - ``(1.0 - exp(U)) / U``

    In addition to the singularity itself, this functions returns a lower and upper bounds as follows:
    - sp is the singularity point: the value of V for which U == 0
    - Vmin is the value of V for which U - U_offset == 0
    - Vmax is the value of V for which U + U_offset == 0:
    :return: (Vmin, Vmax, sp)
    """
    # Create "wildcards", that act as catch-all (*) when matching expressions
    # https://docs.sympy.org/latest/modules/core.html#sympy.core.basic.Basic.match
    P = Wild('P', real=True, exclude=[V])
    Z = Wild('Z', real=True)
    U_wildcard = Wild('U_wildcard', real=True, include=[V])
    SP_wildcard = Wild('SP_wildcard', real=True)

    def check_top_match(m, sp):
        assert m is None or m[Z] != 0
        return m is not None and \
            (sp == m[SP_wildcard] or (isinstance(sp, Float) and 
              isinstance(m[SP_wildcard], Float) and isclose(m[SP_wildcard], sp)))

    # the denominator is all args where a **-1
    numerator = [a for a in expr.args if not _is_negative_power(a)]
    denominator = [Pow(a.args[0], - a.args[1]) for a in expr.args if _is_negative_power(a)]
    
    # Not a division or does not have exp
    if len(denominator) == 0 or len(numerator) == 0 or not expr.has(exp_function):
        return None, None, None
        
    numerator += [Mul(*numerator)]  # pattern could match entire numerator as well as any of its parts
    denominator += [Mul(*denominator)]  # pattern could match entire denominator as well as any of its parts

    # U might be on top, try numerator / denominator and denominator / numerator
    for num, denom in ((numerator, denominator), (denominator, numerator)):
        found_on_top = False
        for d in denom:  # Check arguments in denominator (or numerator)
            (Vmin, Vmax, sp) = None, None, None
            if not d.has(exp_function):
                continue
            # look for exp(U) * -Z + 1.0
            # Note: the -Z match works for finding negative nmbers because the numbers are still in Qauntities
            # Preventing the matching engine from 
            find_U = d.match(exp_function(U_wildcard) * -Z + 1.0)  
            if not find_U:
                find_U = d.match(exp_function(U_wildcard) * Z - 1.0)  # look for exp(U) * Z - 1.0
            if find_U and find_U[Z] > 0:
                # We found a match, since exp(U) * Z == exp(U + log(Z)) we can bring Z into the u expression
                u = (find_U[U_wildcard] + log(find_U[Z]))
                # Find the singularity point by solving u for V==0, excluding irrational results
                sp = solve(u, V)#_solve_tuple(u,V)
                if sp and len(sp) > 0:
                    assert len(sp) == 1, 'the pattern matching should bring numbers added /removed outside the exp, so we should only get 1 solution ' +str(sp)
                    sp = sp[0]
                    Vmin = solve(u - U_offset, V)
                    Vmax = solve(u + U_offset, V)
                    if len(Vmin) > 0 and len (Vmax) > 0:
                        Vmin, Vmax = min(Vmin), max(Vmax)
                    else:  # fallback we can't solve u +/ 1e-7
                        Vmin = sp - U_offset
                        Vmax = sp + U_offset
                    for n in num:  # Check arguments in numerator (or denominator)
                        match = n.match(P * u)  # search for multiple of U
                        found_on_top = match is not None and P in match and match[P] != 0
                        if not found_on_top:
                            match = n.match(Z * V - Z * SP_wildcard)  # search for a multiple of V - sp                        
                            found_on_top = check_top_match(match, sp)
                            if not found_on_top:
                                # search for a exp(multiple of V - sp)
                                match = n.match(exp_function(Z * V - Z * SP_wildcard))
                                found_on_top = check_top_match(match, sp)
                            if found_on_top:  # We've found a match stop looking in the other numerator arguments
                                break
                    if found_on_top:  # found singularity, no need to try further
                        break
        if Vmin is not None and found_on_top:  # found singularity, no need to try further
            break
        else:
            (Vmin, Vmax, sp) = None, None, None

    # Put dummies back in and return the singularity point and range boundries
    return (_float_dummies(Vmin), _float_dummies(Vmax), _float_dummies(sp))


def _fix_expr_parts(expr, V, U_offset, exp_function):
    """
    Removes fixable singularities and replaces them with piecewise expressions.

    :return: either (Vmin, Vmax, singularity point, expression, True)
                    if we have identified a singularity needs to be made but it can't be done yet.
                    For eample if 2 singularities with the same singularity point are added up,
                    a single singularity is constructed instead.
             or ``(None, None, None, fixed expr, fixed_expr has piecewise)`` otherwise

    see :meth:`remove_fixable_singularities for more details` "
    """

    if not expr.has(exp_function):  # Expressions without exp don't have GHK-like equations
        return (None, None, None, expr, False)

    if isinstance(expr, Mul):  # 1 * A --> A (remove unneeded 1 *)
        expr = Mul(*[a for a in expr.args if not str(a) in ('1.0', '1')])

    # Turn Quantity dummies into numbers, to enable analysis
    subs_dict = {d: d.evalf(FLOAT_PRECISION) for d in expr.atoms(Quantity)}
    check_U_expr = expr.xreplace(subs_dict)

    if isinstance(expr, Add):  # A + B + ..
        # The expression is an addition, find singularities in each argument
        new_expr_parts = []
        for a in expr.args:
            new_expr_parts.append(_fix_expr_parts(a, V, U_offset, exp_function))

        # If all arguments have the same singularity point, return 1 singularity with the widest range
        range = [item for (Vmin, Vmax, _, _, _) in new_expr_parts for item in (Vmin, Vmax)]
        if len(new_expr_parts) > 1 and len(set([str(sp) for (_, _, sp, _, _) in new_expr_parts])) == 1 \
                and all(isinstance(b, Quantity) for b in range):
            sp = new_expr_parts[0][2]
            Vmin, Vmax = min(range, key=lambda v: float(str(v))), max(range, key=lambda v: float(str(v)))
            return (Vmin, Vmax, sp, expr, True)
        else:  # Singularity points differ, so create each piecewise seperately and add them up
            expr_parts = []
            is_piecewise = False
            for Vmin, Vmax, sp, ex, has_piecewise in new_expr_parts:
                is_piecewise = is_piecewise or has_piecewise or Vmin is not None
                expr_parts.append(_generate_piecewise(ex, V, sp, Vmin, Vmax) if sp is not None else ex)
            return (None, None, None, Add(*expr_parts), is_piecewise)

    elif isinstance(expr, Pow) and expr.args[1] == -1.0 and len(expr.args) == 2:  # 1/A
        # Find singularities in A and adjust result to represent 1 / A
        Vmin, Vmax, sp, ex, has_piecewise = _fix_expr_parts(expr.args[0], V, U_offset, exp_function)
        has_piecewise = has_piecewise or Vmin is not None
        return (None, None, None, ONE / (_generate_piecewise(ex, V, sp, Vmin, Vmax) if sp is not None else ex),
                has_piecewise)

    elif isinstance(expr, Mul):  # A * B * ...
        (Vmin, Vmax, sp) = _get_singularity(check_U_expr, V, U_offset, exp_function)  # Find the singularity point
        if sp is not None:
            return (Vmin, Vmax, sp, expr, True)
        else:  # Couldn't find singularity, try the expression's arguments
            expr_parts = []
            is_piecewise = False
            for sub_ex in expr.args:
                Vmin, Vmax, sp, ex, has_piecewise = _fix_expr_parts(sub_ex, V, U_offset, exp_function)
                has_piecewise = is_piecewise = is_piecewise or has_piecewise or Vmin is not None
                expr_parts.append(_generate_piecewise(ex, V, sp, Vmin, Vmax) if sp is not None else ex)
            return (None, None, None, Mul(*expr_parts), is_piecewise)

    else:  # Different type of equation e.g a number
        return (None, None, None, expr, False)


def _remove_singularities(expr, V, U_offset=1e-7, exp_function=exp):
    """
    Removes suitable singularities and replaces it with a piecewise.
    :param expr: The expression to analyse ().
    :param V: The voltage variable.
    :param U_offset: Determins the offset either side of U (see get_singularity) for which the fix is used.
    :param exp_function: The function representing exp.
    :return: (bool, expr) as follows: (expr has changed, expr with singularities fixes if appropriate).

    see :meth:`remove_fixable_singularities for more details `."
    """
    if not expr.has(exp_function):
        return False, expr
    (Vmin, Vmax, sp, ex, changed) = _fix_expr_parts(expr, V, U_offset, exp_function)
    ex = _generate_piecewise(ex, V, sp, Vmin, Vmax) if sp is not None else ex
    return changed or Vmin is not None, ex


def remove_fixable_singularities(model, V, modifiable_parameters, U_offset=1e-7, exp_function=exp):
    """
    Finds singularities in the GHK-like equations in the model and replaces them with a piecewise.
    :param model: The cellmlmanip model the model to analyse.
    :param V: The voltage variable.
    :param modifiable_parameters: The variables which are modifiable in the model,
            their defining equations are excluded form the analysis.
    :param U_offset: Determins the offset either side of U for which the fix is used.
    :param exp_function: the function representing exp.
    """
    unprocessed_eqs = {}
    # iterate over sorted variables in the model
    for variable in nx.lexicographical_topological_sort(model.graph, key=str):

        # Get equation
        eq = model.graph.nodes[variable]['equation']

        # Skip variables that have no equation or equations defining parameters or where rhs is a Piecewise
        if eq is not None and not isinstance(eq.rhs, Piecewise) and eq.lhs not in modifiable_parameters:
            unprocessed_eqs[eq.lhs] = eq.rhs.xreplace(unprocessed_eqs)  # store partially evaluated version of the rhs
            changed, new_ex = _remove_singularities(unprocessed_eqs[eq.lhs], V, U_offset=U_offset,
                                                    exp_function=exp_function)
            if changed:  # update equation if the rhs has a singularity that can be fixed
                model.remove_equation(eq)
                model.add_equation(Eq(eq.lhs, new_ex))
                unprocessed_eqs.pop(eq.lhs)
