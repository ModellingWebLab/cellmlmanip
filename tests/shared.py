# Module with cellmlmanip tests
import os

import cellmlmanip


OXMETA = "https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#"


def load_model(name, unit_store=None, skip_singularity_fixes=False):
    """ Parses and returns one of the CellML test models """
    if name[-7:] != '.cellml':
        name += '.cellml'
    return cellmlmanip.load_model(
        os.path.join(os.path.dirname(__file__), 'cellml_files', name),
        unit_store=unit_store, skip_singularity_fixes=skip_singularity_fixes)


def check_left_right_units_equal(unit_store, equality):
    """
    Checks whether the LHS and RHS in a ``sympy.Eq`` have the same units.
    :param unit_store: a :class:`UnitStore`.
    :param equality: A ``sympy.Eq``.
    """
    lhs_units = unit_store.evaluate_units(equality.lhs)
    rhs_units = unit_store.evaluate_units(equality.rhs)
    assert unit_store.is_equivalent(rhs_units, lhs_units), 'Units %s %s != %s %s' % (
        unit_store.format(lhs_units), unit_store.format(lhs_units, True),
        unit_store.format(rhs_units), unit_store.format(rhs_units, True),
    )
