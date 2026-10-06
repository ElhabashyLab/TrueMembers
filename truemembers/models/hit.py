import uuid
from dataclasses import dataclass, field

from truemembers.models.detection_method import DetectionMethod


@dataclass(frozen=True)
class Hit:
    uniprot_id: str
    region_start: int
    region_end: int
    detection_method: DetectionMethod
    unique_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    bitscore: float | None = None
    tm_score_to_target: float | None = None
    ali_coverage: float | None = None
    domains: dict[str, tuple[int, int]] = field(default_factory=dict)
    is_all_member: bool = False
    is_strict_member: bool = False
    first_seen_iteration: int = 0
    structure_path: str | None = None
    sequence_path: str | None = None

    def __post_init__(self) -> None:
        if not self.uniprot_id:
            raise ValueError("Hit.uniprot_id must be a non-empty string")
        if self.region_start < 1:
            raise ValueError(f"Hit.region_start must be >= 1, got {self.region_start}")
        if self.region_end < self.region_start:
            raise ValueError(
                f"Hit.region_end ({self.region_end}) must be >= region_start ({self.region_start})"
            )

    @property
    def region(self) -> tuple[int, int]:
        return self.region_start, self.region_end
