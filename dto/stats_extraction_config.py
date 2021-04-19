from dataclasses import dataclass


@dataclass
class StatsExtractionConfig:
    data_path: str
    output_path: str
    properties: list
    filter: dict
    columns: list
    mode: str = "orig_data"
    token_distribution: bool = True
    plotting: bool = True
    standardisation_config: dict = None
    save_standardised: bool = False
    save_cut_precomputed: bool = False
    token_atom_ratio: bool = True
    count_decorations: bool = False

