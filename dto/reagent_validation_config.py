from dataclasses import dataclass

@dataclass
class ReagentValidationConfig:
    input_smarts: str
    output_smarts: str

