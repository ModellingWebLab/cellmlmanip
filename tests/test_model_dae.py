import os
import pytest

from cellmlmanip import parser


class TestModelDAE:
    """
    Tests if unsupported algebraic models are detected.
    """

    def test_dae(self):

        # Parsing should be OK
        p = parser.Parser(os.path.join(os.path.dirname(__file__), 'cellml_files', '4.algebraic_ode_model.cellml'))
        model = p.parse()

        # But equation graph will raise error
        with pytest.raises(RuntimeError, match='DAEs are not supported'):
            model.get_equation_graph()

