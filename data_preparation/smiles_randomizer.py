import os
from collections import OrderedDict
from typing import List, Dict, Tuple

from reinvent_chemistry import Conversions
from reinvent_chemistry.library_design.attachment_points import AttachmentPoints
from reinvent_chemistry.library_design import BondMaker
from data_preparation.enums.dataframe_columns_enum import DataframeColumnsEnum

import utils.spark as us
from dto import RandomizationConfig


class SmilesRandomizer:
    def __init__(self, configuration: RandomizationConfig):
        self._configuration = configuration
        self._bond_maker = BondMaker()
        self._conversions = Conversions()
        self._columns = DataframeColumnsEnum
        self._attachments = AttachmentPoints()

    def re_label(self, _labels: List[int], _decorations: OrderedDict) -> OrderedDict:
        _reordered_decorations = OrderedDict()
        for _i, v in enumerate(_labels):
            _reordered_decorations[_i] = self.randomize_molecule(_decorations[v])
        return _reordered_decorations

    def randomize_molecule(self, _scaffold_smi: str) -> str:
        _molecule = self._conversions.smile_to_mol(_scaffold_smi)
        result = self._conversions.mol_to_random_smiles(_molecule)
        return result

    def validate_randomization(self, unlabeled: str, relabeled: str, orig: str, decs: str, scaffs: str):
        labeled_result = self._attachments.add_attachment_point_numbers(unlabeled, canonicalize=False)
        molecule = self._bond_maker.join_scaffolds_and_decorations(labeled_result, relabeled)
        reproduction = self._conversions.mol_to_smiles(molecule)
        original = self._conversions.mol_to_smiles(self._conversions.smile_to_mol(orig))
        if reproduction == original:
            pass
        else:
            print(f"{reproduction}, {original} {relabeled} {decs} {scaffs}")

    def row_transformation(self, row: str) -> Dict:

        scaffold_smi, decorations, original = row.split("\t")
        numbered_scaffold = self._attachments.add_attachment_point_numbers(scaffold_smi, canonicalize=False)
        randomized_scaffold = self.randomize_molecule(numbered_scaffold)

        reordered_decorations = OrderedDict()
        decoration_smis = decorations.split("|")

        for i, decoration in enumerate(decoration_smis):
            reordered_decorations[i] = decoration

        labels = self._attachments.get_attachment_points(randomized_scaffold)
        unlabeled_result = self._attachments.remove_attachment_point_numbers(randomized_scaffold)

        relabeled = self.re_label(labels, reordered_decorations)
        relabeled_decorations = '|'.join(relabeled.values())

        if self._configuration.validate_randomization:
            self.validate_randomization(unlabeled=unlabeled_result, relabeled=relabeled_decorations, orig=original,
                                        decs=decorations, scaffs=scaffold_smi)

        return {self._columns.SCAFFOLDS: unlabeled_result, self._columns.DECORATIONS: relabeled_decorations,
                self._columns.ORIGINAL: original}

    def _format_output(self, row: Tuple) -> Tuple:
        return row[self._columns.SCAFFOLDS], row[self._columns.DECORATIONS], row[self._columns.ORIGINAL]

    def run(self):
        SPARK, CONTEXT = us.SparkSessionSingleton.get("create_randomized_smiles")
        os.makedirs(self._configuration.output_path, exist_ok=True)

        for i in range(self._configuration.number_of_files):
            sliced_mols_rdd = CONTEXT.textFile(self._configuration.input_file) \
                .repartition(self._configuration.number_of_partitions) \
                .map(lambda x: self.row_transformation(x)) \
                .persist()

            with open("{}/{:03d}.smi".format(self._configuration.output_path, i), "w+") as out_file:
                for scaff_smi, dec_smi, original in sliced_mols_rdd.map(self._format_output).collect():
                    out_file.write("{}\t{}\t{}\n".format(scaff_smi, dec_smi, original))

