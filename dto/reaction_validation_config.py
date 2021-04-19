from dataclasses import dataclass

@dataclass
class ReactionValidationConfig:
    input_path: str
    output_path: str
    output_smiles_file: str
    reactions_file: str
    failures_limit: int = 10
    num_partitions: int = 1000
    max_cuts: int = 4