import pandas as pd
from reinvent_chemistry import Conversions
from reinvent_chemistry.library_design import BondMaker, AttachmentPoints
from data_preparation.enums.dataframe_columns_enum import DataframeColumnsEnum
from dto.reagent_validation_config import ReagentValidationConfig


class ReagentValidator:
    def __init__(self, config: ReagentValidationConfig):
        self._configuration = config
        self._attachment_points = AttachmentPoints()
        self.chemistry = Conversions()
        self._bond_maker = BondMaker()
        self._columns = DataframeColumnsEnum

    def join_fragments(self, row: pd.DataFrame) -> str:
        scaffold = self._attachment_points.add_attachment_point_numbers(row[self._columns.SCAFFOLDS],
                                                                        canonicalize=False)
        molecule = self._bond_maker.join_scaffolds_and_decorations(scaffold_smi=scaffold,
                                                                   decorations_smi=row[self._columns.DECORATIONS])
        complete_smile = self.chemistry.mol_to_smiles(molecule)
        return complete_smile

    def _load_data(self, path: str) -> pd.DataFrame:
        dataframe = pd.read_csv(path, names=[self._columns.SCAFFOLDS, self._columns.DECORATIONS,
                                             self._columns.ORIGINAL],
                                sep='\t')
        return dataframe

    def _confirm_identity(self, row: pd.DataFrame) -> bool:
        return row['reconstructed'] == row[self._columns.ORIGINAL]

    def run(self):
        data = self._load_data(self._configuration.input_smarts)
        data["reconstructed"] = data.apply(self.join_fragments, axis=1)
        data["identity"] = data.apply(self._confirm_identity, axis=1)
        data.to_csv(self._configuration.output_smarts)
