import unittest
from pyspark.sql import SparkSession
import pandas as pd
from junk.chemical_properties import DatasetProfiling

# not ideal - this requires pyspark

class TestPropertyComputation(unittest.TestCase):
    def setUp(self):
        self.testing_data_single = {'scaffolds': ['[*]C#Cc1ccc2oc(C(C)=O)cc2c1'], 'decorations' : ['*C1(O)CN2CCC1CC2']}
        self.dataset_profiling = DatasetProfiling()

    def _create_testing_pyspark_session(self):
        return SparkSession.builder.master('local[2]').appName('local_testing_pyspark_context').\
            getOrCreate()

    def _get_spark_df(self, data_dict):
        spark = self._create_testing_pyspark_session()
        df = pd.DataFrame(data_dict, index=[1])
        return spark.createDataFrame(df)

    def test_compute_stats(self):
        tested_function = self.dataset_profiling.compute_stats
        df = self._get_spark_df(self.testing_data_single)
        df_processed = tested_function(df, 'scaffolds')
        print(df_processed.show())
        self.assertEqual(1, df_processed.count())

    def test_token_distribution(self):
        tested_function = self.dataset_profiling.token_distribution
        df = self._get_spark_df(self.testing_data_single)
        df_processed = tested_function(df, 'scaffolds')
        print(df_processed)


    # TODO Attachment points etc.
    # def test_number_attachments(self):
    #
    #     self.assertEqual(True, False)
