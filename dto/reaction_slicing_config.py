from dataclasses import dataclass


@dataclass
class ReactionBasedSlicingConfig:
    input_file: str
    output_path: str
    output_smiles_file: str
    conditions_file: str
    reactions_file: str
    max_cuts: int = 4
    number_of_partitions: int = 1000
    validate_randomization: bool = True
