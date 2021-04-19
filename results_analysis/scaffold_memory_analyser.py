from typing import List, Dict, Any
from collections import OrderedDict

from rdkit import Chem
from reinvent_chemistry.library_design import AttachmentPoints

from dto.scaffold_memory_analysis_config import ScaffoldMemoryAnalysisConfig
import pandas as pd
import os

from results_analysis.enums import ScaffoldMemoryPropertyEnum


class ScaffoldMemoryAnalyser:
    def __init__(self, configuration: ScaffoldMemoryAnalysisConfig):
        self.scaffold_memory_paths = configuration.input_paths
        self.output_path = configuration.output_path
        self.top_n_compounds = configuration.n_highest
        self.decoration_analysis_folder = os.path.join(self.output_path, "decoration_analysis")
        self._attachments = AttachmentPoints()
        self._property_enum = ScaffoldMemoryPropertyEnum()
        self._properties = configuration.properties

    def load_scaffold_memory(self, path: str) -> pd.DataFrame:
        try:
            memory_df = pd.read_csv(os.path.join(path, "scaffold_memory.csv"))
            memory_df = memory_df.sort_values("total_score")
            memory_df['run_path'] = path
            memory_df["CanonicalScaffold"] = memory_df.SMILES.apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x)))
            return memory_df

        except FileNotFoundError:
            print(f"Did not find scaffold memory for run {path}")

    def get_average_scores(self, all_memories_df: pd.DataFrame) -> None:
        df_by_run_score = all_memories_df.groupby('run_path')['total_score'].mean()
        df_by_run_score.to_csv(os.path.join(self.output_path, "mean_scaffold_memory_scores.csv"))

    def analyse_whole_compound_overlap(self, all_memories_df: pd.DataFrame, list_of_memory_dfs: List) -> None:
        count_df = all_memories_df.groupby("CanonicalScaffold")['run_path'].count().sort_values(ascending=False)
        overlap = [smiles for smiles in count_df[count_df > 1].index]

        with open(os.path.join(self.output_path, f"overlap_analysis_{self.top_n_compounds}.csv"), "a+") as f:
            heading = "SMILES"
            for run in self.scaffold_memory_paths:
                heading += f"\t{run}"
            f.write(f"{heading}\n")

            for smiles in overlap:
                f.write(f"{smiles}")
                for df in list_of_memory_dfs:
                    try:
                        line = df[df.CanonicalScaffold == smiles]
                        f.write(f'\t{line.Step.item()}')
                    except:
                        f.write("\t NotFound")
                f.write("\n")

    def re_label(self, _labels: List[int], _decorations: OrderedDict) -> OrderedDict:
        _reordered_decorations = OrderedDict()
        for _i, v in enumerate(_labels):
            _reordered_decorations[_i] = self.canonicalise_molecule(_decorations[v])
        return _reordered_decorations

    def canonicalise_molecule(self, _scaffold_smi: str) -> str:
        _molecule = Chem.MolFromSmiles(_scaffold_smi)
        result = Chem.MolToSmiles(_molecule)
        return result

    def row_transformation(self, row: str) -> Dict:
        scaffold_smi = row.split("|")[0]
        decoration_smis = [d for d in row.split("|")[1:]]
        numbered_scaffold = self._attachments.add_attachment_point_numbers(scaffold_smi, canonicalize=False)
        randomized_scaffold = self.canonicalise_molecule(numbered_scaffold)

        reordered_decorations = OrderedDict()

        for i, decoration in enumerate(decoration_smis):
            reordered_decorations[i] = decoration
        labels = self._attachments.get_attachment_points(randomized_scaffold)
        unlabeled_result = self._attachments.remove_attachment_point_numbers(randomized_scaffold)
        relabeled = self.re_label(labels, reordered_decorations)
        relabeled_decorations = '|'.join(relabeled.values())

        return {"CanonicalScaffolds": unlabeled_result, "CanonicalDecorations": relabeled_decorations}

    def decoration_frequency_single(self, scaffold_memory_df: pd.DataFrame) -> Dict:
        decorations_df = pd.DataFrame(scaffold_memory_df.Scaffold.apply(lambda r: self.row_transformation(r))).\
            Scaffold.apply(pd.Series)
        decorations_only = pd.DataFrame(decorations_df.CanonicalDecorations.apply(lambda x: x.split("|")).tolist(),
                                        index=decorations_df.index)
        frequencies_per_attachment = {}
        for col in decorations_only.columns:
            frequencies = decorations_only[col].value_counts().sort_values(ascending=False)
            frequencies_per_attachment.update({str(col): frequencies})

        return frequencies_per_attachment

    def fill_in_decoration_frequency_dict(self, scaffold_memory: pd.DataFrame, memory_id: int,
                                          current_dec_freq_dict: Dict) -> Dict:
        frequency_per_attachment = self.decoration_frequency_single(scaffold_memory)

        for attachment in frequency_per_attachment.keys():
            if attachment in current_dec_freq_dict.keys():
                current_dec_freq_dict[attachment].update({memory_id: frequency_per_attachment[attachment]})
            else:
                current_dec_freq_dict.update({attachment: {memory_id: frequency_per_attachment[attachment]}})
        return current_dec_freq_dict

    def decoration_frequency_heatmap(self, frequency_dict: Dict, attachment_number: Any) -> None:
        decorations_df = pd.DataFrame.from_dict(frequency_dict).fillna(0)
        decorations_df['number_of_occurrences'] = decorations_df.sum(axis=1)
        decorations_df.to_csv(os.path.join(self.decoration_analysis_folder,
                              f"decoration_heatmap{attachment_number}.csv"))

    def run(self) -> None:
        all_scaffold_memories = pd.DataFrame()
        all_df_list = []
        decoration_frequency_dict = {}

        for i, path in enumerate(self.scaffold_memory_paths):
            df = self.load_scaffold_memory(path)
            if self.top_n_compounds > 0:
                df = df.head(self.top_n_compounds)
            all_scaffold_memories = all_scaffold_memories.append(df)
            all_df_list.append(df)

            if self._property_enum.DECORATION_STATS in self._properties:
                print("Updating dictionaries of frequencies")
                if not os.path.exists(self.decoration_analysis_folder):
                    os.mkdir(self.decoration_analysis_folder)
                decoration_frequency_dict = self.fill_in_decoration_frequency_dict(scaffold_memory=df, memory_id=i,
                                                                                   current_dec_freq_dict=
                                                                                   decoration_frequency_dict)

        if self._property_enum.DECORATION_STATS in self._properties:
            print("Producing heatmap")
            for attachment in decoration_frequency_dict.keys():
                self.decoration_frequency_heatmap(decoration_frequency_dict[attachment], attachment)

        if self._property_enum.SCORES in self._properties:
            print("Computing scores per run")
            self.get_average_scores(all_scaffold_memories)

        if self._property_enum.OVERLAP in self._properties:
            print("Computing overlap dataframe")
            self.analyse_whole_compound_overlap(all_scaffold_memories, all_df_list)
