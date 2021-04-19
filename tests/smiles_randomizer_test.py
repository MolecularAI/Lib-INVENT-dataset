import unittest
from collections import OrderedDict

from rdkit import Chem

from data_preparation.smiles_randomizer import SmilesRandomizer
from dto import RandomizationConfig
from data_preparation.enums.dataframe_columns_enum import DataframeColumnsEnum
from tests.fixtures.reaction_fixtures import RELABELING_DICT, SCAFFOLD, ROW


class TestSmilesRandomizer(unittest.TestCase):
    def setUp(self):

        self.relabeling_case = RELABELING_DICT
        # self.scaffolds = SCAFFOLD
        # self.row = ROW
        # self.scaffolds = "*C(=O)C(*)(C)CC(=O)OC1CCC2(*)C(CCC3(C)C2CCC2C4C(*)CCC4(C(=O)OCc4ccccc4)CCC23C)C1(C)C"
        # self.row = "*C(=O)C(*)(C)CC(=O)OC1CCC2(*)C(CCC3(C)C2CCC2C4C(*)CCC4(C(=O)OCc4ccccc4)CCC23C)C1(C)C	*O|*C|*C(C)C|*C	CC(C)C1CCC2(C(=O)OCc3ccccc3)CCC3(C)C(CCC4C5(C)CCC(OC(=O)CC(C)(C)C(=O)O)C(C)(C)C5CCC43C)C12"
        self.row = "[*]C#CC(c1ccc(OCc2cccc(CN([*])[*])c2)cc1)[*]	*C|*C(C)C|*Cc1ccsc1|*CC(=O)O	CC#CC(CC(=O)O)c1ccc(OCc2cccc(CN(Cc3ccsc3)C(C)C)c2)cc1"
        self.scaffolds = "[*]C#CC(c1ccc(OCc2cccc(CN([*])[*])c2)cc1)[*]"
        config = RandomizationConfig("","")
        self.smiles_randomizer = SmilesRandomizer(config)
        self._columns = DataframeColumnsEnum()

    def test_relabel(self):
        ordered_decs = OrderedDict(sorted(self.relabeling_case.items(), key=lambda t: t[1]))
        labels = [v for v in self.relabeling_case.keys()]

        relabeled = self.smiles_randomizer.re_label(_labels=labels, _decorations=ordered_decs)
        self.assertNotEqual(relabeled, ordered_decs)

        for key in labels:
            old = ordered_decs[key]
            new = relabeled[key]
            self.assertEqual(Chem.MolToSmiles(Chem.MolFromSmiles(old)),Chem.MolToSmiles(Chem.MolFromSmiles(new)))

    def test_randomize_molecule(self):
        original = self.scaffolds
        randomized = self.smiles_randomizer.randomize_molecule(original)
        self.assertNotEqual(original, randomized)
        self.assertEqual(Chem.MolToSmiles(Chem.MolFromSmiles(original)),
                         Chem.MolToSmiles(Chem.MolFromSmiles(randomized)))


    def test_row_transformation(self):
        transformed = self.smiles_randomizer.row_transformation(self.row)
        orig_scaffold, orig_decorations, orig_full = self.row.split("\t")
        self.assertNotEqual(orig_scaffold, transformed[self._columns.SCAFFOLDS])
        self.assertNotEqual(orig_decorations, transformed[self._columns.DECORATIONS])
        self.assertEqual(Chem.MolToSmiles(Chem.MolFromSmiles(orig_full)),
                         Chem.MolToSmiles(Chem.MolFromSmiles(transformed[self._columns.ORIGINAL])))
        self.assertEqual(Chem.MolToSmiles(Chem.MolFromSmiles(orig_scaffold)),
                         Chem.MolToSmiles(Chem.MolFromSmiles(transformed[self._columns.SCAFFOLDS])))
