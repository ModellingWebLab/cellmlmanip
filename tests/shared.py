# Module with cellmlmanip tests
import logging


logger = logging.getLogger(__name__)


def check_cmeta_ids(model):
    """Checks that every variable in a model with a cmeta_id is a source variable."""
    is_okay = True
    for variable in model.graph:
        if variable.is_Derivative:
            variable = variable.free_symbols.pop()
        if variable.cmeta_id:
            if variable != variable.assigned_to:
                is_okay = False
                logger.critical('%s has cmeta id but is assigned to %s',
                                variable.dummy, variable.assigned_to)
    return is_okay


def check_dummy_assignment(model):
    """Every variable in the model should be assigned to itself or another variabe.
    The source variable must be assigned to itself."""
    is_okay = True
    for variable in model.graph:
        if variable.is_Derivative:
            variable = variable.free_symbols.pop()
        # either the variable is assigned to itself
        if variable == variable.assigned_to:
            continue

        # or the variable is assigned to a source variable
        source = variable.assigned_to

        # the source dummy must be assigned to itself
        if source.assigned_to != source:
            is_okay = False
            logger.critical('%s is assigned to %s, which is assigned to %s',
                            variable,
                            source,
                            source.assigned_to)
    return is_okay
