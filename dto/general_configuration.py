from dataclasses import dataclass


@dataclass
class GeneralConfiguration:
    run_type: str
    parameters: dict
