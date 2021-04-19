import os
import re
from typing import List, Dict

import matplotlib.pyplot as plt
import pandas as pd
import pyspark.sql as ps
import pyspark.sql.functions as psf
import pyspark.sql.types as pst
from pyspark.sql import SparkSession
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Descriptors import ExactMolWt
from reinvent_chemistry.enums import FilterTypesEnum

from reinvent_chemistry.standardization.rdkit_standardizer import RDKitStandardizer
from reinvent_chemistry.standardization.filter_configuration import FilterConfiguration

import utils.spark as us
from data_preparation.enums import DataframeColumnsEnum, StatsExtractionEnum
from data_preparation.enums.purging_enum import PurgingEnum
from dto.stats_extraction_config import StatsExtractionConfig


class StatsExtractor:
    """
    A large class used jointly for data preprocessing and analysis, since the methods are typically shared.
    The data has to either be a single column of full compounds or a sliced dataset containing scaffolds, decorations
    and original compounds. The default mode is evaluation of the single column data.
    The computed properties can be saved for analysis or used as a basis for filtering the dataset.
    The available properties are: ['mol_wts', 'num_rings', 'num_aromatic_rings', 'num_atoms', 'hbond_donors', 'hbond_acceptors']
    """
    def __init__(self, configuration: StatsExtractionConfig):
        self._filters = FilterTypesEnum

        self._columns = DataframeColumnsEnum
        self._stats = StatsExtractionEnum
        self._purging = PurgingEnum
        self._configuration = configuration
        standardisation_config_dict = self._configuration.standardisation_config
        standardisation_config = [FilterConfiguration(name=name, parameters=params)
                                  for name, params in standardisation_config_dict.items()]

        dec_separator = self._stats.DECORATION_SEPARATOR_TOKEN
        attachment_token = self._stats.ATTACHMENT_POINT_TOKEN
        self._mol_wts_udf = psf.udf(lambda x: ExactMolWt(Chem.MolFromSmiles(x)), pst.FloatType())
        self._num_rings_udf = psf.udf(lambda x: rdMolDescriptors.CalcNumRings(Chem.MolFromSmiles(x)),
                                      pst.IntegerType())
        self._num_atoms_udf = psf.udf(lambda x: Chem.MolFromSmiles(x).GetNumHeavyAtoms(), pst.IntegerType())
        self._num_aromatic_rings_udf = psf.udf(lambda x: rdMolDescriptors.CalcNumAromaticRings(Chem.MolFromSmiles(x)),
                                               pst.IntegerType())
        self._hbond_donors_udf = psf.udf(lambda x: rdMolDescriptors.CalcNumHBD(Chem.MolFromSmiles(x)),
                                         pst.IntegerType())
        self._hbond_acceptors_udf = psf.udf(lambda x: rdMolDescriptors.CalcNumHBA(Chem.MolFromSmiles(x)),
                                            pst.IntegerType())
        self._hetero_atom_ratio_udf = psf.udf(lambda x: len([atom for atom in Chem.MolFromSmiles(x).GetAtoms()
                                                            if atom.GetAtomicNum() == 6])/Chem.MolFromSmiles(x).
                                              GetNumHeavyAtoms(),
                                              pst.FloatType())
        self._make_canonical_udf = psf.udf(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x)), pst.StringType())
        self._standardise_smiles_udf = psf.udf(lambda x: RDKitStandardizer(standardisation_config, None).apply_filter(x),
                                               pst.StringType())
        pattern = self._stats.REGEX_TOKENS
        self.regex = re.compile(pattern)
        self._tokeniser_udf = psf.udf(self.regex.findall, pst.ArrayType(pst.StringType()))
        self._decoration_split_udf = psf.udf(lambda x: x.split(dec_separator),
                                             pst.ArrayType(pst.StringType()))
        self._count_decorations_udf = psf.udf(lambda s: list(s).count(attachment_token),
                                              pst.IntegerType())

    def _number_attachment_points(self, df: ps.DataFrame) -> pd.DataFrame:
        attachments_histogram = df.withColumn(self._stats.NUMBER_OF_DECORATIONS,
                                              self._count_decorations_udf(self._columns.SCAFFOLDS)).\
            groupBy(self._stats.NUMBER_OF_DECORATIONS).count().toPandas().\
            set_index(self._stats.NUMBER_OF_DECORATIONS).sort_index()
        return attachments_histogram

    def _count_decorations(self, df: ps.DataFrame) -> None:
        print('Working on decoration count')
        df = self._number_attachment_points(df)
        save_folder = os.path.join(self._configuration.output_path, 'data')

        save_path = os.path.join(save_folder, f'{self._stats.NUMBER_OF_DECORATIONS}.csv')
        df.to_csv(save_path)

        if self._configuration.plotting:
            self._plot_and_save(df, self._stats.NUMBER_OF_DECORATIONS, bar=True)

    def smiles_tokenizer(self, smi: str) -> List[str]:
        tokens = [token for token in self.regex.findall(smi)]
        if smi != ''.join(tokens):
            print(smi, ''.join(tokens))
            raise AssertionError
        return tokens

    def _compute_single_property(self, df: ps.DataFrame, property: str, column) -> ps.DataFrame:
        f_name = f'_{property}_udf'
        fn = getattr(self, f_name)
        df = df.withColumn(property, fn(column))
        return df

    def _compute_all_stats(self, df: ps.DataFrame, methods: List[str], column) -> ps.DataFrame:
        for method in methods:
            df = self._compute_single_property(df, method, column)
        return df

    def update_token_dictionary(self, token_dictionary: Dict, character_list: List[str]) -> Dict:
        for token in character_list:
            try:
                token_dictionary[token] += 1
            except KeyError:
                token_dictionary[token] = 1
        return token_dictionary

    def token_distribution(self, df: ps.DataFrame, column: str) -> pd.DataFrame:
        token_dict = dict()
        for row in df.rdd.toLocalIterator():
            row_list = self.smiles_tokenizer(row[column])
            token_dict = self.update_token_dictionary(token_dict, row_list)
        out_df = pd.DataFrame.from_dict(token_dict, orient='index').rename(columns={0: 'counts'})\
            .sort_values('counts', ascending=False)
        return out_df

    def _plot_and_save(self, df: pd.DataFrame, name: str, bar: bool = False) -> None:
        if bar:
            df.plot.bar()
        else:
            df.plot()
        save_path_image = os.path.join(self._configuration.output_path, f'plots/{name}.png')
        plt.savefig(save_path_image)

    def _make_folders(self):
        save_folder = os.path.join(self._configuration.output_path, 'data')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        if self._configuration.plotting:
            save_folder = os.path.join(self._configuration.output_path, 'plots')
            if not os.path.exists(save_folder):
                os.makedirs(save_folder)

    def _compute_stats_for_property(self, df: ps.DataFrame, property: str, column: str) -> None:
        print(f'Working on {property}')
        pandas_df = df.groupBy(property).agg(psf.count(property).alias('counts')). \
            toPandas().set_index(property).sort_index()
        save_path = os.path.join(self._configuration.output_path, f'data/{property}_{column}.csv')
        pandas_df.to_csv(save_path)
        if self._configuration.plotting:
            self._plot_and_save(pandas_df, f"{property}_{column}", bar=False)

    def _compute_token_distribution_stats(self, df: ps.DataFrame, column: str) -> None:
        print('Working on token distribution')
        token_df = self.token_distribution(df, column)
        save_path = os.path.join(self._configuration.output_path, f'data/{self._stats.TOKEN_DISTRIBUTION}_{column}.csv')
        token_df.to_csv(save_path)
        if self._configuration.plotting:
            self._plot_and_save(token_df, name=f"{self._stats.TOKEN_DISTRIBUTION}_{column}", bar=True)

    def _compute_token_atom_ratio(self, df: ps.DataFrame, column: str) -> None:
        token_numbers = df.select(column)
        token_numbers = token_numbers.withColumn("token_list", self._tokeniser_udf(column))
        token_numbers = token_numbers.withColumn(self._stats.NUMBER_OF_TOKENS, psf.size("token_list"))

        self._compute_stats_for_property(token_numbers, self._stats.NUMBER_OF_TOKENS, column)
        if self._stats.NUMBER_OF_ATOMS not in self._configuration.properties:
            print("Missing number of atoms. Can't get ratios.")
        else:
            df = df.select(column, self._stats.NUMBER_OF_ATOMS)
            df = df.join(token_numbers, on=[column])
            df = df.withColumn(self._stats.TOKEN_ATOM_RATIO,
                               df[self._stats.NUMBER_OF_TOKENS]/df[self._stats.NUMBER_OF_ATOMS])
            self._compute_stats_for_property(df, self._stats.TOKEN_ATOM_RATIO, column)

    def _filter_data(self, df: ps.DataFrame) -> ps.DataFrame:
        for prop, condition in self._configuration.filter.items():
            if prop not in df.columns:
                print(f"cannot filter on {prop} if it is not computed.")
                pass
            else:
                rule = " ".join([prop] + [self._purging.MAPPING[condition[0]]] + [str(condition[1])])
                print(rule)
                df = df.where(rule)
        return df

    def prepare_single_column_data(self) -> ps.DataFrame:
        self._make_folders()
        SPARK, CONTEXT = us.SparkSessionSingleton.get("preprocessor")
        data = SPARK.read.csv(self._configuration.data_path)
        data.show()
        data = data.select(self._standardise_smiles_udf("_c0").alias(self._columns.ORIGINAL)).dropDuplicates() \
            .where(f"{self._columns.ORIGINAL} is not null")
        data.show()
        if self._configuration.save_standardised:
            data.coalesce(1).write.format("text").option("header", "false").mode("append"). \
                save(self._configuration.output_path)
        return data

    def prepare_sliced_data(self, data: ps.DataFrame,  column: str) -> ps.DataFrame:
        self._make_folders()
        if column == self._columns.DECORATIONS:
            data = data.select(psf.explode(self._decoration_split_udf(self._columns.DECORATIONS)).
                               alias(self._columns.DECORATIONS))

        elif column == self._columns.SCAFFOLDS:
            if self._configuration.count_decorations:
                self._count_decorations(df=data)
        return data

    def run_computations(self, data: ps.DataFrame, column: str) -> None:
        data_computed = self._compute_all_stats(data, self._configuration.properties, column)
        data_computed = self._filter_data(data_computed)

        if self._configuration.save_cut_precomputed:
            print("****about to save****")
            data_computed.select(column).coalesce(1).write.format("text").option("header", "false"). \
                mode("append").save(self._configuration.output_path)

        print("****about to compute stats****")
        for _method in self._configuration.properties:
            self._compute_stats_for_property(data_computed, _method, column=column)

        if self._configuration.token_distribution:
            self._compute_token_distribution_stats(data_computed, column)

        if self._configuration.token_atom_ratio:
            self._compute_token_atom_ratio(data_computed, column)

    def run(self):

        if self._configuration.mode == self._stats.SINGLE_COLUMN_DATA_MODE:
            data = self.prepare_single_column_data()
            self.run_computations(data, self._columns.ORIGINAL)

        elif self._configuration.mode == self._stats.SLICED_DATA_MODE:
            SPARK = SparkSession.builder.appName('data_stats').getOrCreate()
            data = SPARK.read.options(delimiter='\t').csv(self._configuration.data_path). \
                withColumnRenamed('_c0', self._columns.SCAFFOLDS).withColumnRenamed('_c1', self._columns.DECORATIONS). \
                withColumnRenamed('_c2', self._columns.ORIGINAL)

            for column in self._configuration.columns:
                print(f"*****working on column: {column}******")
                data = self.prepare_sliced_data(data, column)
                self.run_computations(data, column)
        else:
            print("Unknown preprocessing mode.")
