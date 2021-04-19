from dataclasses import dataclass

@dataclass
class DecorationSimilarityAnalyserConfig:
    decoration_frequency_df_path: str
    output_path: str
    standardise: bool = True