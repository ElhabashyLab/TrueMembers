from enum import Enum


class StopReason(str, Enum):
    CONVERGED = "CONVERGED"
    MAX_ITERATIONS = "MAX_ITERATIONS"
    SEEDS_EXHAUSTED = "SEEDS_EXHAUSTED"
    FAILED = "FAILED"
