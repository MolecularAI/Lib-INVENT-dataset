from dataclasses import dataclass
from typing import List


@dataclass
class ScaffoldMemoryAnalysisConfig:
    input_paths: List
    output_path: str
    properties: List
    n_highest: int = -1