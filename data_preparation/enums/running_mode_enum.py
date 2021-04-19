

class RunningModeEnum:
    RANDOMIZATION = "data_randomization"
    REACTION_BASED_SLICING = "reaction_based_slicing"
    DUPLICATE_REMOVAL = "duplicate_removal"
    STATS_EXTRACTION = "stats_extraction"
    FILE_SHUFFLING = "file_shuffling"
    REACTION_VALIDATION = "reaction_validation"
    REAGENT_VALIDATION = "reagent_validation"
    VALIDATION_SCAFFOLD_SELECTION = "validation_scaffold_selection"
    TENSORBOARD_LOG_EXTRACTION = "tensorboard_log_extraction"
    SCAFFOLD_MEMORY_ANALYSIS = "scaffold_memory_analysis"
    DECORATION_SIMILARITY = "decoration_similarity"

    VALIDATION_SET_SIMILARITIES = "scaffold_similarity"
    VALIDATION_SET_FILTERING = "filtering"
    VALIDATION_SET_SLICED = "sliced_split"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
