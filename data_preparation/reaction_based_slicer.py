import json
from typing import List

import pyspark.sql as ps
import pyspark.sql.functions as psf
from reinvent_chemistry.library_design import FragmentReactions, AttachmentPoints
from reinvent_chemistry.library_design.dtos import SlicingConditionsDTO, FilteringConditionDTO
from reinvent_chemistry.library_design.fragment_reaction_slice_enumerator import FragmentReactionSliceEnumerator
from data_preparation.enums.dataframe_columns_enum import DataframeColumnsEnum


import utils.chem as uc
import utils.spark as us
from dto import ReactionBasedSlicingConfig


class ReactionBasedSlicer:

    def __init__(self, configuration: ReactionBasedSlicingConfig):
        self._configuration = configuration
        self._columns = DataframeColumnsEnum

    def _to_smiles_rows(self, row: ps.Row) -> str:
        return "{}\t{}\t{}".format(row[self._columns.SCAFFOLDS], row[self._columns.DECORATIONS],
                                   row[self._columns.ORIGINAL])

    @staticmethod
    def enumerate(row: ps.Row, enumerator: FragmentReactionSliceEnumerator, max_cuts: int) -> List[ps.Row]:
        attachments = AttachmentPoints()
        fields = row.split("\t")
        smiles = fields[0]
        mol = uc.to_mol(smiles)
        out_rows = []
        if mol:
            for sliced_mol in enumerator.enumerate(mol, cuts=max_cuts):
                row_dict = {
                    DataframeColumnsEnum.SCAFFOLDS:
                        attachments.remove_attachment_point_numbers(sliced_mol.scaffold_smiles),
                    DataframeColumnsEnum.DECORATIONS: sliced_mol.decorations_smiles,
                    DataframeColumnsEnum.ORIGINAL: sliced_mol.original_smiles,
                    DataframeColumnsEnum.MAX_CUTS: max_cuts}
                out_rows.append(ps.Row(**row_dict))
        return out_rows

    def slice_db(self, enumerator: FragmentReactionSliceEnumerator) -> ps.DataFrame:
        spark, sc = us.SparkSessionSingleton.get("slice_db")
        enumeration_df = spark.createDataFrame(
            sc.textFile(self._configuration.input_file)
                .repartition(self._configuration.number_of_partitions)
                .flatMap(lambda f: self.enumerate(f, enumerator=enumerator, max_cuts=self._configuration.max_cuts))) \
            .groupBy(self._columns.SCAFFOLDS, self._columns.DECORATIONS) \
            .agg(psf.first(self._columns.MAX_CUTS).alias(self._columns.MAX_CUTS),
                 psf.first(self._columns.ORIGINAL).alias(self._columns.ORIGINAL)) \
            .persist()

        if self._configuration.output_path:
            enumeration_df.write.parquet(self._configuration.output_path)
        return enumeration_df

    def run(self):
        scaffold_conditions = []
        decoration_conditions = []
        if self._configuration.conditions_file:
            with open(self._configuration.conditions_file, "r") as json_file:
                data = json.load(json_file)
                conditions = SlicingConditionsDTO(**data)
                scaffold_conditions = [FilteringConditionDTO(**condition) for condition in conditions.scaffold]
                decoration_conditions = [FilteringConditionDTO(**condition) for condition in conditions.decoration]

        chemistry = FragmentReactions()

        if self._configuration.reactions_file:
            with open(self._configuration.reactions_file, "r") as reactions_file:
                lines = reactions_file.read().splitlines()
        reactions = chemistry.create_reactions_from_smirks(lines)

        _enumerator = FragmentReactionSliceEnumerator(reactions, scaffold_conditions, decoration_conditions)
        slice_df = self.slice_db(_enumerator)

        if self._configuration.output_smiles_file:
            with open(self._configuration.output_smiles_file, "w+") as smiles_file:
                for row in slice_df.rdd.map(self._to_smiles_rows).toLocalIterator():
                    smiles_file.write("{}\n".format(row))

