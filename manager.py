from dacite import from_dict
from results_analysis.scaffold_memory_analyser import ScaffoldMemoryAnalyser
from data_preparation.stats_extractor import StatsExtractor
from data_preparation.reaction_validator import ReactionValidator
from data_preparation.reagent_validator import ReagentValidator
from data_preparation.smiles_randomizer import SmilesRandomizer
from data_preparation.duplicate_finder import DuplicateFinder
from data_preparation.reaction_based_slicer import ReactionBasedSlicer
from data_preparation.file_shuffler import FileShuffler
from results_analysis.tensorboard_logs import TensorboardStatsExtractor

from dto import GeneralConfiguration, RandomizationConfig, ReactionBasedSlicingConfig, \
    DuplicateRemovalConfig, StatsExtractionConfig, FileShufflingConfig, ReactionValidationConfig
from dto.reagent_validation_config import ReagentValidationConfig
from data_preparation.enums.running_mode_enum import RunningModeEnum
from dto.scaffold_memory_analysis_config import ScaffoldMemoryAnalysisConfig
from dto.tensorboard_extraction_config import TensorboardStatsExtractionConfig


class Manager:
    def __init__(self, configuration: GeneralConfiguration):
        self._configuration = configuration
        self._running_mode = RunningModeEnum()

    def _data_randomization(self):
        config = from_dict(data_class=RandomizationConfig, data=self._configuration.parameters)
        randomizer = SmilesRandomizer(config)
        randomizer.run()

    def _reaction_based_slicing(self):
        config = from_dict(data_class=ReactionBasedSlicingConfig, data=self._configuration.parameters)
        slicer = ReactionBasedSlicer(config)
        slicer.run()

    def _duplicate_removal(self):
        config = from_dict(data_class=DuplicateRemovalConfig, data=self._configuration.parameters)
        remover = DuplicateFinder(config)
        remover.run()

    def _file_shuffling(self):
        config = from_dict(data_class=FileShufflingConfig, data=self._configuration.parameters)
        shuffler = FileShuffler(config)
        shuffler.run()

    def _reaction_validation(self):
        config = from_dict(data_class=ReactionValidationConfig, data=self._configuration.parameters)
        validator = ReactionValidator(config)
        validator.run()

    def _reagent_validation(self):
        config = from_dict(data_class=ReagentValidationConfig, data=self._configuration.parameters)
        validator = ReagentValidator(config)
        validator.run()

    def _stats_extractor(self):
        config = from_dict(data_class=StatsExtractionConfig, data=self._configuration.parameters)
        preprocessor = StatsExtractor(config)
        preprocessor.run()

    def _tensorboard_extractor(self):
        config = from_dict(data_class=TensorboardStatsExtractionConfig, data=self._configuration.parameters)
        tensorboard_extractor = TensorboardStatsExtractor(config)
        tensorboard_extractor.run()

    def _scaffold_memory_analyser(self):
        config = from_dict(data_class=ScaffoldMemoryAnalysisConfig, data=self._configuration.parameters)
        analyser = ScaffoldMemoryAnalyser(config)
        analyser.run()

    def run(self):
        registry = {self._running_mode.RANDOMIZATION: self._data_randomization,
                    self._running_mode.REACTION_BASED_SLICING: self._reaction_based_slicing,
                    self._running_mode.DUPLICATE_REMOVAL: self._duplicate_removal,
                    self._running_mode.STATS_EXTRACTION: self._stats_extractor,
                    self._running_mode.FILE_SHUFFLING: self._file_shuffling,
                    self._running_mode.REACTION_VALIDATION: self._reaction_validation,
                    self._running_mode.REAGENT_VALIDATION: self._reagent_validation,
                    self._running_mode.TENSORBOARD_LOG_EXTRACTION: self._tensorboard_extractor,
                    self._running_mode.SCAFFOLD_MEMORY_ANALYSIS: self._scaffold_memory_analyser}

        run_type = registry.get(self._configuration.run_type, lambda: TypeError)
        run_type()
