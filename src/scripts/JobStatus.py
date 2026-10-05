import enum

class JobStatus(enum.Enum):
    PENDING = "PENDING" # Job has been submitted but not yet started
    RUNNING = "RUNNING" # Job is currently running
    COMPLETED = "COMPLETED" # Job has finished successfully
    HARVESTED = "HARVESTED" # Job has completed or failed but got processed and results have been harvested
    FAILED = "FAILED" # Job has finished with an error
    