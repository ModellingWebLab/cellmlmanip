import os

from cellmlmanip import parser


class TestParser(object):
    @staticmethod
    def get_test_cellml_filepath():
        return os.path.join(os.path.dirname(__file__), "cellml_files", "test_simple_odes.cellml")

    def test_loader(self):
        p = parser.Parser(TestParser.get_test_cellml_filepath())
        model = p.parse()
        print(model)

