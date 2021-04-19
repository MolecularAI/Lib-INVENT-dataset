import json
import unittest

from reinvent_chemistry.library_design import FragmentReactions, FragmentReactionSliceEnumerator
from reinvent_chemistry.library_design.dtos import FilteringConditionDTO

from data_preparation.enums import DataframeColumnsEnum
from data_preparation.reaction_based_slicer import ReactionBasedSlicer
from dto import ReactionBasedSlicingConfig
from tests.fixtures.reaction_fixtures import REACTION_SMIRKS_ALL, COMPOUND


class TestReactionBasedSlicing(unittest.TestCase):
    def setUp(self):
        config = ReactionBasedSlicingConfig("","","","","")
        self.scaffold_conditions = [FilteringConditionDTO(name='molecular_weight', min=0)]
        self.decoration_conditions = [FilteringConditionDTO(name='molecular_weight', min=0)]
        self.smirks = REACTION_SMIRKS_ALL
        self.compound = COMPOUND
        self.slicer = ReactionBasedSlicer(config)
        self.chemistry = FragmentReactions()
        self._columns = DataframeColumnsEnum

    def _get_enumerator(self):
        chemistry = FragmentReactions()
        reactions = chemistry.create_reactions_from_smirks(self.smirks)
        enumerator = FragmentReactionSliceEnumerator(reactions, self.scaffold_conditions, self.decoration_conditions)
        return enumerator

    def test_enumerate(self):
        enumerate_fn = self.slicer.enumerate
        enumerator = self._get_enumerator()
        out = enumerate_fn(self.compound, enumerator, max_cuts=4)[0].asDict()
        self.assertEqual(out[self._columns.ORIGINAL], self.compound)
        self.assertEqual(out[self._columns.SCAFFOLDS].count('*'), len(out[self._columns.DECORATIONS].split('|')))
        self.assertGreaterEqual(out['max_cuts'], out[self._columns.SCAFFOLDS].count('*'))

