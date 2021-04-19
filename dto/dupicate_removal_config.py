from dataclasses import dataclass


@dataclass
class DuplicateRemovalConfig:
    data_path: str
    output_path: str
    save: bool = True
    remove: bool = True