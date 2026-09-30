from __future__ import annotations

from dataclasses import dataclass

from src.Hit import Hit, DetectionMethod


DEFAULT_TMSCORE_THRESHOLD = 0.6


@dataclass(frozen=True)
class Member(Hit):
    tm_score_threshold: float = DEFAULT_TMSCORE_THRESHOLD

    def __post_init__(self):
        if self.qtmscore < self.tm_score_threshold:
            raise ValueError(
                f"Cannot create Member: qtmscore {self.qtmscore} is below "
                f"threshold {self.tm_score_threshold}"
            )

    @classmethod
    def from_hit(cls, hit: Hit, tm_score_threshold: float = DEFAULT_TMSCORE_THRESHOLD) -> Member:
        if hit.qtmscore < tm_score_threshold:
            raise ValueError(
                f"Hit qtmscore {hit.qtmscore} is below threshold {tm_score_threshold}"
            )
        return cls(
            unique_id=hit.unique_id,
            uniprot_id=hit.uniprot_id,
            region_start=hit.region_start,
            region_end=hit.region_end,
            sequence_path=hit.sequence_path,
            structure_path=hit.structure_path,
            detection_method=hit.detection_method,
            qtmscore=hit.qtmscore,
            domains=list(hit.domains),
            tm_score_threshold=tm_score_threshold,
        )