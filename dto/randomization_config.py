from dataclasses import dataclass


@dataclass
class RandomizationConfig:
    input_file: str
    output_path: str
    number_of_files: int = 1
    number_of_partitions: int = 1000
    validate_randomization: bool = True
