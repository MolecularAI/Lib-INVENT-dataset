import os
from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
import pandas as pd
from dto.tensorboard_extraction_config import TensorboardStatsExtractionConfig


class TensorboardStatsExtractor:
    def __init__(self, configuration: TensorboardStatsExtractionConfig):
        self.configuration = configuration

    def run(self):
        df_list = []
        for directory in self.configuration.log_folders:
            list_log_folders = os.listdir(directory)
            dir_name = os.path.basename(os.path.dirname(directory))
            list_log_folders = [d for d in list_log_folders if os.path.isdir(os.path.join(directory, d))]

            for tb_output_folder in list_log_folders:
                print("working on:", tb_output_folder)
                x = EventAccumulator(path=os.path.join(directory, tb_output_folder))
                x.Reload()
                x.FirstEventTimestamp()
                keys = self.configuration.properties
                print_out_dict = {}

                for i in range(len(keys)):
                    print_out_dict.update({f"{dir_name}_{keys[i]}_{tb_output_folder}":
                                               [e.value for e in x.Scalars(keys[i])]})

                df = pd.DataFrame(data=print_out_dict)
                df_list.append(df)

        complete_df = pd.concat(df_list, axis=1)
        complete_df.to_csv(os.path.join(self.configuration.output_path, self.configuration.output_name))
