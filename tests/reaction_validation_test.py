import unittest

from reinvent_chemistry.library_design.failing_reactions_enumerator import FailingReactionsEnumerator

from data_preparation.reaction_validator import ReactionValidator
from dto import ReactionValidationConfig
from reinvent_chemistry.library_design import FragmentReactions

from tests.fixtures.reaction_fixtures import WORKING_REACTION, FAILED_REACTION_ROW, FAILED_REACTION


class ReactionValidatorTest(unittest.TestCase):
    def setUp(self):
        config = ReactionValidationConfig("","","","")
        self.validator = ReactionValidator(config)

        self.okay_row = WORKING_REACTION
        self.failed_row = FAILED_REACTION_ROW
        self.second_failed_reaction = FAILED_REACTION


    def _get_enumerator(self, lines):
        chemistry = FragmentReactions()
        reactions = chemistry.create_reactions_from_smirks(lines)
        enumerator = FailingReactionsEnumerator(reactions)
        return enumerator

    def test_collect_failures_single_row(self):
        ''' Test to check that a failed reaction is detected. Need the original compound too because that serves as an input
        for the collect_failures function.'''
        tested = self.validator.collect_failures
        failed_reaction, failed_compound = self.failed_row.split()
        lines = [failed_reaction, self.okay_row]
        enumerator = self._get_enumerator(lines)
        result = tested(failed_compound, enumerator)
        self.assertEqual(result[0][0], failed_compound)
        self.assertEqual(result[0][1], failed_reaction)
        self.assertEqual(len(result), 1)

    def test_collect_failures_one_compound_multiple_reactions(self):
        ''' Test to check multiple failures of the same slicing are detected.'''
        tested = self.validator.collect_failures
        failed_reaction, failed_compound = self.failed_row.split()
        lines = [failed_reaction, self.okay_row, self.second_failed_reaction]
        enumerator = self._get_enumerator(lines)
        result = tested(failed_compound, enumerator)
        self.assertEqual(result[0][0], failed_compound)
        self.assertEqual(result[1][0], failed_compound)
        self.assertEqual(set([result[0][1], result[1][1]]), set([failed_reaction, self.second_failed_reaction]))
        self.assertEqual(len(result), 2)
