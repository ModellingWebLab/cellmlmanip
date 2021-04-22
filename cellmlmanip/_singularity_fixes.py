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
    Mul,
    Piecewise,
    Pow,
    S,
    Wild,
    exp,
    log,
    solveset,
)
from sympy.codegen.rewriting import ReplaceOptim, optimize
from sympy.sets import Intersection

from .model import FLOAT_PRECISION, Quantity


# For P^n make sure n is passed as int if it is actually a whole number
_POW_OPT = ReplaceOptim(lambda p: p.is_Pow and (isinstance(p.exp, Float) or isinstance(p.exp, float))
                        and float(p.exp).is_integer(),
                        lambda p: Pow(p.base, int(float(p.exp))))


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
    """Replaces any sympy ``Float`` objects in ``expr`` with (dimensionless) :class:`cellmlmanip.model.Quantity` objects."""
    return expr.xreplace({f: Quantity(f, 'dimensionless') for f in expr.atoms(Float)})


def _is_negative_power(expr):
    """
    Checks if ``expr`` is a sympy ``Pow`` with a negative exponent.
    
    If the exponent contains variables the sign of the exponent can't be determined, and ``False`` is returned.
    """
    try:
        return isinstance(expr, Pow) and bool(expr.args[1] < 0)
    except TypeError:  # if this is a power with variables still in it we can't determine if it's negative
        pass
    return False


def _solve_real(u, V):
    """Gives the values of V for which u == 0 """
    u = optimize(u, (_POW_OPT, ))  # make sure powers of ints are represented as ints
    result = solveset(u, V, domain=S.Reals)
    # The resul could be an intersection with Reals, if a custom exp function is used. We assume the result is real.
    if isinstance(result, Intersection) and result.args[0] == S.Reals:
        result = result.args[1]
    return result


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
    singularities = []

    def check_U_match(m, sp):
        assert m is None or m[Z] != 0
        return m is not None and \
            (sp == m[SP_wildcard] or (isinstance(sp, Float) and
             isinstance(m[SP_wildcard], Float) and isclose(m[SP_wildcard], sp)))

    # find all fractions and separate numerator and denominator
    # the denominator is all args where a **-1
    numerator = [a for a in expr.args if not _is_negative_power(a)]
    denominator = [Pow(a.args[0], - a.args[1]) for a in expr.args if _is_negative_power(a)]

    # Second check: must be a fraction
    if len(denominator) == 0 or len(numerator) == 0 or not expr.has(exp_function):
        return []

    numerator += [Mul(*numerator)]  # pattern could match entire numerator as well as any of its parts
    denominator += [Mul(*denominator)]  # pattern could match entire denominator as well as any of its parts

    # U might be on top, try numerator / denominator and denominator / numerator
    for fraction_part_1, fraction_part_2 in ((numerator, denominator), (denominator, numerator)):
        found_on_top = False
        for fp2 in fraction_part_2:  # Check arguments in denominator (or numerator)
            (Vmin, Vmax, sp) = None, None, None
            if not fp2.has(exp_function):
                continue
            # look for exp(U) * -Z + 1.0
            # Works as we later require Z to be positive
            find_U = fp2.match(exp_function(U_wildcard) * -Z + 1.0)
            if not find_U:
                find_U = fp2.match(exp_function(U_wildcard) * Z - 1.0)  # look for exp(U) * Z - 1.0

            # Z should be positive, we replace free symbols by 1 to be able to evaluate the sign of Z
            if find_U and find_U[Z].xreplace({s: 1.0 for s in find_U[Z].free_symbols}) > 0:
                # We found a match for exp(U) * Z -1 or exp(U) * -Z +1,
                # since exp(U) * Z == exp(U + log(Z)) we can bring Z into the u expression
                # Note: the top check (check_U_match) will not require Z to be positive (just not 0)
                u = (find_U[U_wildcard] + log(find_U[Z]))

                # Find the singularity point by solving u for V==0, excluding irrational results
                singularity_points = _solve_real(u, V)

                # Find Vmin, Vmax for the fix range
                Vmin_points = _solve_real(u - U_offset, V)
                Vmax_points = _solve_real(u + U_offset, V)
                for sp in singularity_points:
                    if len(Vmin_points) == len(Vmax_points) == 1:
                        Vmin, Vmax = tuple(Vmin_points)[0], tuple(Vmax_points)[0]
                    else:  # multiple solutions for Vmin / Vmax we don't know which to use, so revert to fixed range
                        Vmin, Vmax = sp - U_offset, sp + U_offset

                    for fp1 in fraction_part_1:  # Check arguments in numerator (or denominator)
                        match = fp1.match(P * u)  # search for multiple of U
                        found_on_top = match is not None and P in match and match[P] != 0
                        if not found_on_top:
                            match = fp1.match(Z * V - Z * SP_wildcard)  # search for a multiple of V - sp
                            found_on_top = check_U_match(match, sp)
                            if not found_on_top:
                                # search for a exp(multiple of V - sp)
                                match = fp1.match(exp_function(Z * V - Z * SP_wildcard))
                                found_on_top = check_U_match(match, sp)
                        if found_on_top:  # We've found a match stop looking in the other numerator arguments
                            break
                    if found_on_top:  # found singularity
                        # If this isngularity was previously found, adjust Vmin/Vmax to widest range
                        for sing in singularities:
                            if sing[2] == sp:
                                # if Vmin & Vmax are the same it's the same singularity
                                if Vmin == sing[0] and Vmax == sing[1]:
                                    break
                                elif all(len(v.free_symbols) == 0 for v in (Vmin, Vmax, sing[0], sing[1])):
                                    # if Vmin/Vmax have no variables we can pick min/max
                                    sing[0] = min(sing[0], sing[1], Vmin, Vmax)
                                    sing[1] = max(sing[0], sing[1], Vmin, Vmax)
                                    break
                        else:  # if it's a new singularity add to the list
                            singularities.append([Vmin, Vmax, sp])
        if not found_on_top:
            (Vmin, Vmax, sp) = None, None, None

    # Put dummies back in and return the singularity point and range boundries
    return [(_float_dummies(Vmin), _float_dummies(Vmax), _float_dummies(sp)) for (Vmin, Vmax, sp) in singularities]


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
        singularities = _get_singularity(check_U_expr, V, U_offset, exp_function)  # Find the singularity point
        range = [item for (Vmin, Vmax, _) in singularities for item in (Vmin, Vmax)]
        # If all arguments have the same singularity point, return 1 singularity with the widest range
        if len(singularities) == 1:
            Vmin, Vmax, sp = singularities[0]
            return (Vmin, Vmax, sp, expr, True)
        elif len(singularities) > 0:  # Singularity points differ, so generate a nested piecewise
            for Vmin, Vmax, sp in singularities[1:]:
                expr = _generate_piecewise(expr, V, sp, Vmin, Vmax)
            Vmin, Vmax, sp = singularities[0]
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
