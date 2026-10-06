from dataclasses import dataclass

from truemembers.models.job_status import JobStatus


@dataclass
class JobSummary:
    job_id: str
    iteration: int
    seed_id: str
    seed_uniprot_id: str
    seed_region_start: int
    seed_region_end: int
    status: JobStatus = JobStatus.COMPLETED
    attempts: int = 1
    sequence_hits: int = 0
    structure_hits: int = 0
    deduped_hits: int = 0
    downloads: int = 0
    structural_pass: int = 0
    new_strict_members: int = 0
    duration_s: float = 0.0
    error: str | None = None
