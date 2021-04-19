import unittest

from data_preparation.stats_extractor import StatsExtractor
from dto import StatsExtractionConfig

from data_preparation.enums.stats_extraction_enum import StatsExtractionEnum
from tests.fixtures.stats_extraction_fixtures import SCAFFOLD1, TOKEN_DICT_SCAFFOLD1

# TODO this needs more detailed testing.
class StatsExtractorTest(unittest.TestCase):
    def setUp(self):
        self.scaffold = SCAFFOLD1
        config = StatsExtractionConfig
        self.exceptions = StatsExtractionEnum.TOKEN_EXCEPTIONS
        self.extractor = StatsExtractor(config)

    def test_token_count(self):
        count = self.extractor.token_count(self.scaffold)
        self.assertEqual(count, 34)

    def test_update_token_dictionary(self):
        orig_dict = dict()
        char_list = list(self.scaffold)
        updated_dict = self.extractor.update_token_dictionary(orig_dict, char_list)
        expected_dict = TOKEN_DICT_SCAFFOLD1
        self.assertEqual(updated_dict, expected_dict)


