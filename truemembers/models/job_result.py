from dataclasses import dataclass, field

from truemembers.models.hit import Hit


@dataclass
class JobResult:
    seed: Hit
    hits: list[Hit] = field(default_factory=list)
    sequence_hits: int = 0
    structure_hits: int = 0
    deduped_hits: int = 0
    downloads: int = 0
    structural_pass: int = 0
    error: str | None = None

    @property
    def failed(self) -> bool:
        return self.error is not None
