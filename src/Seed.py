from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from src.Hit import Hit


@dataclass(frozen=True)
class Seed:
    hit: Hit
    job_config_path: Path
    slurm_script_path: Path

    @property
    def unique_id(self) -> str:
        return self.hit.unique_id

    @property
    def region(self) -> tuple[int, int]:
        return self.hit.get_target_region

    @property
    def qtmscore(self) -> float:
        return self.hit.qtmscore

    @property
    def is_member(self) -> bool:
        return self.hit.is_member