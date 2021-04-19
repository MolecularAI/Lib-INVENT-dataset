import os

import pyspark.sql as ps
import pyspark.sql.functions as psf
from pyspark.sql import SparkSession

from data_preparation.enums.dataframe_columns_enum import DataframeColumnsEnum
from dto import DuplicateRemovalConfig


class DuplicateFinder:
    def __init__(self, configuration: DuplicateRemovalConfig):
        self._configuration = configuration
        self._columns = DataframeColumnsEnum

    def _to_smiles_rows(self, row: ps.Row) -> str:
        return "{}\t{}\t{}".format(row[self._columns.SCAFFOLDS], row[self._columns.DECORATIONS],
                                   row[self._columns.ORIGINAL])

    def _find_duplicates(self, df: ps.DataFrame, config: DuplicateRemovalConfig) -> None:
        df = df.groupBy(self._columns.SCAFFOLDS, self._columns.DECORATIONS).agg(psf.count(self._columns.ORIGINAL).
                                                                                alias('count'),
                                                                                psf.collect_list(
                                                                                    self._columns.ORIGINAL).
                                                                                alias(self._columns.ORIGINAL))
        df = df.filter(df['count'] > 1)
        number_repetitions = df.count()
        print(f'The number of duplicates in the dataset is {number_repetitions}')

        if config.save:
            destination = os.path.join(config.output_path, 'duplicates.smi')
            with open(destination, "w+") as smiles_file:
                for row in df.rdd.map(self._to_smiles_rows).toLocalIterator():
                    smiles_file.write("{}\n".format(row))

    def run(self):
        spark = SparkSession.builder.appName('data_stats').getOrCreate()
        data = spark.read.options(delimiter='\t').csv(self._configuration.data_path). \
            withColumnRenamed('_c0', DataframeColumnsEnum.SCAFFOLDS). \
            withColumnRenamed('_c1', DataframeColumnsEnum.DECORATIONS). \
            withColumnRenamed('_c2', DataframeColumnsEnum.ORIGINAL)

        self._find_duplicates(data, config=self._configuration)

        if self._configuration.remove:
            removed = data.dropDuplicates()

            destination = os.path.join(self._configuration.output_path, 'no_duplicate.smi')
            with open(destination, "w+") as smiles_file:
                for row in removed.rdd.map(self._to_smiles_rows).toLocalIterator():
                    smiles_file.write("{}\n".format(row))
