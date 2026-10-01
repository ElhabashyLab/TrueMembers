from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from src.Hit import Hit


@dataclass(frozen=True)
class Target:
    hit: Hit
    roi_structure_path: Path

    @property
    def unique_id(self) -> str:
        return self.hit.unique_id

    @property
    def uniprot_id(self) -> str:
        return self.hit.uniprot_id

    @property
    def region(self) -> tuple[int, int]:
        return self.hit.get_target_region