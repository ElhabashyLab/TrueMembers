from dataclasses import dataclass, field
from typing import Any

from truemembers.models.duration import duration_seconds
from truemembers.models.iteration_summary import IterationSummary
from truemembers.models.job_summary import JobSummary
from truemembers.models.stop_reason import StopReason


@dataclass
class RunResult:
    target_uid: str
    target_region_start: int
    target_region_end: int
    stop_reason: StopReason
    run_dir: str
    iterations: list[IterationSummary] = field(default_factory=list)
    jobs: list[JobSummary] = field(default_factory=list)
    total_all_members: int = 0
    total_strict_members: int = 0
    started_at: str | None = None
    finished_at: str | None = None
    config_snapshot: dict[str, Any] = field(default_factory=dict)

    @property
    def duration_s(self) -> float | None:
        return duration_seconds(self.started_at, self.finished_at)
