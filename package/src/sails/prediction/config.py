import dataclasses
import enum


@dataclasses.dataclass
class Configuration:
    """Configuration for Sails"""

    use_gpu: bool = False
    n_threads: int | None = None
    disable_progress_bar: bool = True
    compute_entire_unit_cell: bool = True
    compute_variance: bool = False
    use_raw_values: bool = False
    spacing: float = 0.7
    box_size: int = 128
    channels: int = 2
    overlap: int = 64


class MapType(enum.Enum):
    """Map types for sails, i.e. model will output 1 for protein..."""

    glycan: int = 1
    protein: int = 2
