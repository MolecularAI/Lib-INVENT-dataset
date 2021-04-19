from typing import List

import pyspark.sql as ps
from reinvent_chemistry.library_design import FragmentReactions
from reinvent_chemistry.library_design.failing_reactions_enumerator import FailingReactionsEnumerator
from data_preparation.enums.dataframe_columns_enum import DataframeColumnsEnum

import utils.chem as uc
import utils.spark as us
from dto import ReactionValidationConfig


class ReactionValidator:
    def __init__(self, config: ReactionValidationConfig):
        self._columns = DataframeColumnsEnum
        self.configuration = config

    def _to_smiles_rows(self, row: ps.Row) -> str:
        return "{}\t{}".format(row[self._columns.REACTION], row[self._columns.ORIGINAL])

    def collect_failures(self, row: ps.Row, enumerator: FailingReactionsEnumerator) -> List[ps.Row]:
        fields = row.split("\t")
        smiles = fields[0]
        mol = uc.to_mol(smiles)
        out_rows = []
        if mol:
            for failed_reaction in enumerator.enumerate(mol, failures_limit=self.configuration.failures_limit):
                row_dict = {
                    self._columns.REACTION: failed_reaction.reaction_smirks,
                    self._columns.ORIGINAL: failed_reaction.molecule_smiles
                }
                print("found failed reaction")
                out_rows.append(ps.Row(**row_dict))
                if self.configuration.failures_limit <= len(out_rows):
                    break
        return out_rows

    def get_failed_reaction_df(self) -> ps.DataFrame:
        SPARK, SC = us.SparkSessionSingleton.get("slice_db")
        chemistry = FragmentReactions()

        if self.configuration.reactions_file:
            with open(self.configuration.reactions_file, "r") as reactions_file:
                lines = reactions_file.read().splitlines()

        reactions = chemistry.create_reactions_from_smirks(lines)
        enumerator = FailingReactionsEnumerator(reactions)

        return SPARK.createDataFrame(SC.textFile(self.configuration.input_path).repartition(
            self.configuration.num_partitions).flatMap(lambda x: self.collect_failures(x, enumerator))).persist()

    def run(self):
        df = self.get_failed_reaction_df()

        if self.configuration.output_path:
            df.write.parquet(self.configuration.output_path)

        if self.configuration.output_smiles_file:

            with open(self.configuration.output_smiles_file, "w+") as smiles_file:
                for row in df.rdd.map(self._to_smiles_rows).toLocalIterator():
                    smiles_file.write("{}\n".format(row))

