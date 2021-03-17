from math import isclose

import networkx as nx
from cellmlmanip.model import FLOAT_PRECISION, Quantity
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
    log,
    solveset,
    exp
)


ONE = Quantity(1.0, 'dimensionless')


def _generate_piecewise(vs, ve, sp, ex, V):
    """Generates a piecewsie for expression based on the singularity point sp and vmin (vs) and vmax (ve) """
    def f(Vx, e):
        return e.xreplace({V: Vx})

    if vs is None:  # This shouldn't be apiecewise since we have no Vmin / Vmax
        return ex

    return Piecewise(*[(f(vs, ex) + ((V - vs) / (ve - vs)) * (f(ve, ex) - f(vs, ex)),
                        Abs(V - sp) < Abs((ve - vs) / 2)), (ex, True)])


def _get_U(expr, V, U_offset, exp_function):
    '''Finds U in ghk equations these are of one of the the following forms where U is an expression over V:
       - `U / (exp(U) - 1.0)`
       - `U / (1.0 - exp(U))`
       - `(exp(U) - 1.0) / U`
       - `(1.0 - exp(U)) / U`
       '''
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
        '''Turns flaots back into Quantity dummies to be in line with the rest of the sympy equations'''
        if expr is None:
            return None
        return expr.xreplace({f: Quantity(f, 'dimensionless') for f in expr.atoms(Float)})

    # the denominator is all args where a **-1
    numerator = tuple(a for a in expr.args if not isinstance(a, Pow) or a.args[1] != -1.0)
    denominator = tuple(a.args[0] for a in expr.args if isinstance(a, Pow) and a.args[1] == -1.0)

    # Not a devision or does not have exp
    if len(denominator) == 0 or len(numerator) == 0 or not expr.has(exp_function):
        return None, None, None

    (vs, ve, sp) = None, None, None
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
                            # we found a singularity point, now find vs,ve (or vmin, vmax)
                            # by solving for U_offset either side of U
                            find_v_low = solveset(u + U_offset, V)
                            find_v_up = solveset(u - U_offset, V)
                            find_v_low = tuple(find_v_low)
                            find_v_up = tuple(find_v_up)
                            if find_v_low and find_v_up:
                                assert len(find_v_low) == len(find_v_up) == len(sp) == 1, \
                                    'Expecting exactly 1 solution for singularity point'
                                (vs, ve, sp) = (find_v_low[0], find_v_up[0], sp[0])
                    except TypeError:
                        pass  # Result could be 'ConditionSet' which is not iterable and not Real

                if vs is not None:  # check top
                    for n in num:  # Check arguments in numerator (or denominator)
                        match = n.match(Z * V - Z * SP_wildcard)  # search for a multiple of V - sp
                        found_on_top = check_top_match(match, sp)
                        if not found_on_top:
                            # search for a exp(multiple of V - sp)
                            match = n.match(exp_function(Z * V - Z * SP_wildcard))
                            found_on_top = check_top_match(match, sp)
                            if not found_on_top:
                                # A few equations don't play ball with the SP wildcard
                                # We can find those by adding sp in the pattern directly
                                match = n.match(Z * V - Z * sp)
                                found_on_top = match is not None and Z in match and Z != 0
                                if not found_on_top:
                                    match = n.match(exp_function(Z * V - Z * sp))
                                    found_on_top = match is not None and Z in match and Z != 0
                        if found_on_top:  # We've found a match stop looking in the other numerator arguments
                            break
                    if found_on_top:  # found singularity, no need to try further
                        break
                    else:
                        (vs, ve, sp) = None, None, None
        if vs is not None and found_on_top:  # found singularity, no need to try further
            break
        else:
            (vs, ve, sp) = None, None, None

    # Put dummies back in and return the singularity point and range bundries
    return (float_dummies(vs), float_dummies(ve), float_dummies(sp))


def _new_expr_parts(expr, V, U_offset, exp_function):
    """Removes suitable singularities and replaces it with a piecewise, returning (vs, ve, sp, has_singularity)"""

    if isinstance(expr, Mul):  # 1 * A --> A (remove unneeded 1 *)
        expr = Mul(*[a for a in expr.args if not str(a) in ('1.0', '1')])

    # Turn Quantity dummies into numbers, to enable analysis
    subs_dict = {d: d.evalf(FLOAT_PRECISION) for d in expr.atoms(Quantity)}
    check_U_expr = expr.xreplace(subs_dict)

    if not expr.has(exp_function):  # Expressions without exp don't have GHK equations
        return (None, None, None, expr, False)

    elif isinstance(expr, Add):  # A + B + ..
        # The expression is an addition, find singularities in each argument
        new_expr_parts = []
        for a in expr.args:
            new_expr_parts.append(_new_expr_parts(a, V, U_offset, exp_function))

        # If all arguments have the same singularity point, return 1 singularity with the widest range
        range = [item for (vs, ve, _, _, _) in new_expr_parts for item in (vs, ve)]
        if len(new_expr_parts) > 1 and len(set([str(sp) for (_, _, sp, _, _) in new_expr_parts])) == 1 \
                and all(isinstance(b, Quantity) for b in range):
            sp = new_expr_parts[0][2]
            vs, ve = min(range, key=lambda v: float(str(v))), max(range, key=lambda v: float(str(v)))
            return (vs, ve, sp, expr, True)
        else:  # Singularity points differ, so create each piecewise seperately and add them up
            expr_parts = []
            is_piecewise = False
            for vs, ve, sp, ex, has_piecewise in new_expr_parts:
                is_piecewise = is_piecewise or has_piecewise or vs is not None
                expr_parts.append(_generate_piecewise(vs, ve, sp, ex, V))
            return (None, None, None, Add(*expr_parts), is_piecewise)

    elif isinstance(expr, Pow) and expr.args[1] == -1.0 and len(expr.args) == 2:  # 1/A
        # Find singularities in A and adjust result to represent 1 / A
        vs, ve, sp, ex, has_piecewise = _new_expr_parts(expr.args[0], V, U_offset, exp_function)
        has_piecewise = has_piecewise or vs is not None
        return (None, None, None, ONE / _generate_piecewise(vs, ve, sp, ex, V), has_piecewise)

    elif isinstance(expr, Mul):  # A * B * ...
        (vs, ve, sp) = _get_U(check_U_expr, V, U_offset, exp_function)  # Find the singularity point
        if vs is not None:
            return (vs, ve, sp, expr, True)
        else:  # Couldn't find singularity, try the expression's arguments
            expr_parts = []
            is_piecewise = False
            for sub_ex in expr.args:
                vs, ve, sp, ex, has_piecewise = _new_expr_parts(sub_ex, V, U_offset, exp_function)
                has_piecewise = is_piecewise = is_piecewise or has_piecewise or vs is not None
                expr_parts.append(_generate_piecewise(vs, ve, sp, ex, V))
            return (None, None, None, Mul(*expr_parts), is_piecewise)

    else:  # Different type of equation e.g a number
        return (None, None, None, expr, False)


def new_expr(expr, V, U_offset=1e-7, exp_function=exp):
    """Removes suitable singularities and replaces it with a piecewise.
    :param: expr the expression to analyse ().
    :param: V the volatge variable
    :param: U_offset determins the offset either side of U for which the fix is used
    :param: exp_function the function representing exp
    :param: optimize whether or not to apply optimisations to the new expression
    :return: expr with singularities fixes if appropriate
    """
    if not expr.has(exp_function):
        return False, expr
    (vs, ve, sp, ex, changed) = _new_expr_parts(expr, V, U_offset, exp_function)
    ex = _generate_piecewise(vs, ve, sp, ex, V)
    return changed or vs is not None, ex


def fix_singularity_equations(model, V, modifiable_parameters, U_offset=1e-7, exp_function=exp):
    """Finds singularities in the GHK equations in the model and replaces them with a piecewise.
    :param: the cellmlmanip model the model to analyse.
    :param: V the volatge variable
    :param: modifiable_parameters the variables which are modifiable in the model,
            their defining equations are excluded form the analysis
    :param: U_offset determins the offset either side of U for which the fix is used
    :param: exp_function the function representing exp
    :param: optimize whether or not to apply optimisations to the new expression
    """
    unprocessed_eqs = {}
    # iterate over sorted variables in the model
    for variable in nx.lexicographical_topological_sort(model.graph, key=str):

        # Get equation
        eq = model.graph.nodes[variable]['equation']

        # Skip variables that have no equation or equations defining parameters or where rhs is a Piecewise
        if eq is not None and not isinstance(eq.rhs, Piecewise) and eq.lhs not in modifiable_parameters:
            unprocessed_eqs[eq.lhs] = eq.rhs.xreplace(unprocessed_eqs)  # store partially evaluated version of the rhs
            changed, new_ex = new_expr(unprocessed_eqs[eq.lhs], V, U_offset=U_offset, exp_function=exp_function)
            if changed:  # update equation if the rhs has a singularity that can be fixed
                model.remove_equation(eq)
                model.add_equation(Eq(eq.lhs, new_ex))
                unprocessed_eqs.pop(eq.lhs)
