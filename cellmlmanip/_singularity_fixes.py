"""
``remove_fixable_singularities`` specifies a method to remove fixable singularities from a model
"""
from math import isclose

import networkx as nx
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
)

import cellmlmanip
from cellmlmanip.model import FLOAT_PRECISION


ONE = cellmlmanip.Quantity(1.0, 'dimensionless')


def _generate_piecewise(expr, V, sp, Vmin, Vmax):
    """Generates a new (piecewise) expression based on expr with a linear interpolation when V is between Vmin and Vmax

       The returned expression between Vmin and Vmax is::
       f(vs) + (f(ve) - f(vs)) * (V - vs) / (ve - vs)
       whenever sp is not Null, returns the original expression otherwise.

       :param: The expression to turn into a piecewise
       :param: V the voltage variable.
       :param: the value of the singularity point.
       :param: the value of the lower bound.
       :param: the value of the upper bound.
       :return: expr with a singularity exception around sp, if sp is not None else return the original expr
    """
    if sp is None:  # This shouldn't be apiecewise since we have no Vmin / Vmax
        return expr

    f_Vmin = expr.xreplace({V: Vmin})
    f_Vmax = expr.xreplace({V: Vmax})
    return Piecewise((f_Vmin + ((V - Vmin) / (Vmax - Vmin)) * (f_Vmax - f_Vmin),
                      Abs(V - sp) < Abs((Vmax - Vmin) / 2)), (expr, True))


def _get_singularity(expr, V, U_offset, exp_function):
    """ Finds singularities in equations of the form::
       - ``U / (exp(U) - 1.0)``
       - ``U / (1.0 - exp(U))``
       - ``(exp(U) - 1.0) / U``
       - ``(1.0 - exp(U)) / U``

       In addition to the singularity itself, this functions returns a lower and upper bounds::
       :return: (Vmin, Vmax, singularity point)
       """
    Z = Wild('Z', real=True)
    U_wildcard = Wild('U_wildcard', real=True, include=[V])
    SP_wildcard = Wild('SP_wildcard', real=True)

    def check_bottom_match(m):
        return m is not None and U_wildcard in m and Z in m and Z != 0

    def check_top_match(m, sp):
        return m is not None and Z in m and Z != 0 and (sp == m[SP_wildcard]
                                                        or (isinstance(sp, Float)
                                                            and isinstance(m[SP_wildcard], Float)
                                                            and isclose(m[SP_wildcard], sp)))

    def float_dummies(expr):
        """Turns floats back into Quantity dummies to be in line with the rest of the sympy equations"""
        if expr is None:
            return None
        return expr.xreplace({f: cellmlmanip.Quantity(f, 'dimensionless') for f in expr.atoms(Float)})

    # the denominator is all args where a **-1
    numerator = tuple(a for a in expr.args if not isinstance(a, Pow) or a.args[1] != -1.0)
    denominator = tuple(a.args[0] for a in expr.args if isinstance(a, Pow) and a.args[1] == -1.0)

    # Not a devision or does not have exp
    if len(denominator) == 0 or len(numerator) == 0 or not expr.has(exp_function):
        return None, None, None

    (Vmin, Vmax, sp) = None, None, None
    # U might be on top, try numerator / denominator and denominator / numerator
    for num, denom in ((numerator, denominator), (denominator, numerator)):
        found_on_top = False
        for d in denom:  # Check arguments in denominator (or numerator)
            if d.has(exp_function):
                find_U = d.match(exp_function(U_wildcard) * -Z + 1.0)  # look for exp(U) * -Z + 1.0
                if not check_bottom_match(find_U):
                    find_U = d.match(exp_function(U_wildcard) * Z - 1.0)  # look for exp(U) * Z - 1.0
                if check_bottom_match(find_U):
                    # We found a match, since exp(U) * Z == exp(U + log(Z)) we can bring Z into the u expression
                    u = (find_U[U_wildcard] + log(find_U[Z]))
                    try:
                        # Find the singularity point by solving u for V==0, excluding irratinal results
                        sp = tuple(filter(lambda s: not s.has(I), solveset(u, V)))
                        if sp:
                            # we found a singularity point, now find Vmin, Vmax
                            # by solving for U_offset either side of U
                            find_v_low = solveset(u + U_offset, V)
                            find_v_up = solveset(u - U_offset, V)
                            find_v_low = tuple(find_v_low)
                            find_v_up = tuple(find_v_up)
                            assert find_v_low and find_v_up and len(find_v_low) == len(find_v_up) == len(sp) == 1, \
                                'Expecting exactly 1 singularity point '
                            (Vmin, Vmax, sp) = (find_v_low[0], find_v_up[0], sp[0])
                    except TypeError:
                        pass  # Result could be 'ConditionSet' which is not iterable and not Real

                if Vmin is not None:  # check top
                    for n in num:  # Check arguments in numerator (or denominator)
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
                    else:
                        (Vmin, Vmax, sp) = None, None, None
        if Vmin is not None and found_on_top:  # found singularity, no need to try further
            break
        else:
            (Vmin, Vmax, sp) = None, None, None

    # Put dummies back in and return the singularity point and range bundries
    return (float_dummies(Vmin), float_dummies(Vmax), float_dummies(sp))


def _fixe_expr_parts(expr, V, U_offset, exp_function):
    """Removes fixable singularities and replaces them with piecewise expressions.

    :return: either (Vmin, Vmax, singularity point, expression, True)
                    if we have identified a singularity needs to be made but it can't be done yet.
                    For eample if 2 singularities with the same singularity point are added up,
                    a single singularity is constructed instead.
             or ``(None, None, None, fixed expr, fixed_expr has piecewise)`` otherwise

    see :meth:`fix_singularity_equations for more details` "
    """

    if isinstance(expr, Mul):  # 1 * A --> A (remove unneeded 1 *)
        expr = Mul(*[a for a in expr.args if not str(a) in ('1.0', '1')])

    # Turn Quantity dummies into numbers, to enable analysis
    subs_dict = {d: d.evalf(FLOAT_PRECISION) for d in expr.atoms(cellmlmanip.Quantity)}
    check_U_expr = expr.xreplace(subs_dict)

    if not expr.has(exp_function):  # Expressions without exp don't have GHK-like equations
        return (None, None, None, expr, False)

    elif isinstance(expr, Add):  # A + B + ..
        # The expression is an addition, find singularities in each argument
        new_expr_parts = []
        for a in expr.args:
            new_expr_parts.append(_fixe_expr_parts(a, V, U_offset, exp_function))

        # If all arguments have the same singularity point, return 1 singularity with the widest range
        range = [item for (Vmin, Vmax, _, _, _) in new_expr_parts for item in (Vmin, Vmax)]
        if len(new_expr_parts) > 1 and len(set([str(sp) for (_, _, sp, _, _) in new_expr_parts])) == 1 \
                and all(isinstance(b, cellmlmanip.Quantity) for b in range):
            sp = new_expr_parts[0][2]
            Vmin, Vmax = min(range, key=lambda v: float(str(v))), max(range, key=lambda v: float(str(v)))
            return (Vmin, Vmax, sp, expr, True)
        else:  # Singularity points differ, so create each piecewise seperately and add them up
            expr_parts = []
            is_piecewise = False
            for Vmin, Vmax, sp, ex, has_piecewise in new_expr_parts:
                is_piecewise = is_piecewise or has_piecewise or Vmin is not None
                expr_parts.append(_generate_piecewise(ex, V, sp, Vmin, Vmax))
            return (None, None, None, Add(*expr_parts), is_piecewise)

    elif isinstance(expr, Pow) and expr.args[1] == -1.0 and len(expr.args) == 2:  # 1/A
        # Find singularities in A and adjust result to represent 1 / A
        Vmin, Vmax, sp, ex, has_piecewise = _fixe_expr_parts(expr.args[0], V, U_offset, exp_function)
        has_piecewise = has_piecewise or Vmin is not None
        return (None, None, None, ONE / _generate_piecewise(ex, V, sp, Vmin, Vmax), has_piecewise)

    elif isinstance(expr, Mul):  # A * B * ...
        (Vmin, Vmax, sp) = _get_singularity(check_U_expr, V, U_offset, exp_function)  # Find the singularity point
        if sp is not None:
            return (Vmin, Vmax, sp, expr, True)
        else:  # Couldn't find singularity, try the expression's arguments
            expr_parts = []
            is_piecewise = False
            for sub_ex in expr.args:
                Vmin, Vmax, sp, ex, has_piecewise = _fixe_expr_parts(sub_ex, V, U_offset, exp_function)
                has_piecewise = is_piecewise = is_piecewise or has_piecewise or Vmin is not None
                expr_parts.append(_generate_piecewise(ex, V, sp, Vmin, Vmax))
            return (None, None, None, Mul(*expr_parts), is_piecewise)

    else:  # Different type of equation e.g a number
        return (None, None, None, expr, False)


def _remove_singularities(expr, V, U_offset=1e-7, exp_function=exp):
    """Removes suitable singularities and replaces it with a piecewise.
    :param: expr the expression to analyse ().
    :param: V the voltage variable
    :param: U_offset determins the offset either side of U (see get_singularity) for which the fix is used
    :param: exp_function the function representing exp
    :return: (bool, expr) as follows: (expr has changed, expr with singularities fixes if appropriate)

    see :meth:`fix_singularity_equations for more details `"
    """
    if not expr.has(exp_function):
        return False, expr
    (Vmin, Vmax, sp, ex, changed) = _fixe_expr_parts(expr, V, U_offset, exp_function)
    ex = _generate_piecewise(ex, V, sp, Vmin, Vmax)
    return changed or Vmin is not None, ex


def remove_fixable_singularities(model, V, modifiable_parameters, U_offset=1e-7, exp_function=exp):
    """Finds singularities in the GHK-like equations in the model and replaces them with a piecewise.
    :param: the cellmlmanip model the model to analyse.
    :param: V the voltage variable
    :param: modifiable_parameters the variables which are modifiable in the model,
            their defining equations are excluded form the analysis
    :param: U_offset determins the offset either side of U for which the fix is used
    :param: exp_function the function representing exp
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
