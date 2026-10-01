from __future__ import annotations

import uuid
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


class DetectionMethod(Enum):
    JACKHMMER = "sequence"
    FOLDSEEK = "structure"
    BOTH = "sequence_and_structure"


@dataclass(frozen=True)
class Hit:
    uniprot_id: str
    region_start: int
    region_end: int
    sequence_path: Path
    structure_path: Path
    detection_method: DetectionMethod
    qtmscore: float
    is_member: bool = False
    domains: list[dict[str, tuple[int, int]]] = field(default_factory=list)
    unique_id: str = field(default_factory=lambda: uuid.uuid4().hex)

    def has_domain(self, domain_name: str) -> bool:
        return any(domain_name in d for d in self.domains)

    def domain_regions(self) -> dict[str, tuple[int, int]]:
        regions: dict[str, tuple[int, int]] = {}
        for d in self.domains:
            regions.update(d)
        return regions

    @property
    def get_target_region(self) -> tuple[int, int]:
        return (self.region_start, self.region_end)
    