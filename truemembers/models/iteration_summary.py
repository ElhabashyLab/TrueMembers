from dataclasses import dataclass

from truemembers.models.duration import duration_seconds


@dataclass
class IterationSummary:
    iteration: int
    seeds_submitted: int
    jobs_completed: int
    jobs_failed: int
    new_strict_members: int
    total_strict_members: int
    growth: float = 0.0
    started_at: str | None = None
    finished_at: str | None = None

    @property
    def duration_s(self) -> float | None:
        return duration_seconds(self.started_at, self.finished_at)
