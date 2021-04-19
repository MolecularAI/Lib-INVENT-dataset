from dataclasses import dataclass


@dataclass
class FileShufflingConfig:
    directory: str
    n_permutations: int = 1
    save_permutations: bool = False
