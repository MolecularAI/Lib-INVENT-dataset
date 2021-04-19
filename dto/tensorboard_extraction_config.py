from dataclasses import dataclass
from typing import List


@dataclass
class TensorboardStatsExtractionConfig:
    log_folders: List
    output_path: str
    output_name: str
    properties: List