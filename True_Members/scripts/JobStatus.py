import enum

class JobStatus(enum.Enum):
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    HARVESTED = "HARVESTED"
    FAILED = "FAILED"
    