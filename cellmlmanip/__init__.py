"""Main module for loading, parsing and manipulating CellML models"""
from ._config import __version__, __version_int__, version  # noqa
from .main import load_model  # noqa
from .model import Quantity  # noqa
from ._singularity_fixes import remove_fixable_singularities
