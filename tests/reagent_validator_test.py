import unittest

from rdkit import Chem

from data_preparation.reagent_validator import ReagentValidator
from data_preparation.enums import DataframeColumnsEnum
from dto.reagent_validation_config import ReagentValidationConfig
from tests.fixtures.reaction_fixtures import SCAFFOLD, DECORATIONS, ORIGINAL


class ReagentValidatorTest(unittest.TestCase):
    def setUp(self):
        config = ReagentValidationConfig("","")
        self.validator = ReagentValidator(config)
        self._columns = DataframeColumnsEnum

        self.row = {self._columns.SCAFFOLDS: SCAFFOLD,
                    self._columns.DECORATIONS: DECORATIONS,
                    self._columns.ORIGINAL: ORIGINAL
                    }


    def test_join_fragments(self):
        reconstructed = self.validator.join_fragments(self.row)
        self.assertEqual(Chem.MolToSmiles(Chem.MolFromSmiles(self.row[self._columns.ORIGINAL])),
                         Chem.MolToSmiles(Chem.MolFromSmiles(reconstructed)))
