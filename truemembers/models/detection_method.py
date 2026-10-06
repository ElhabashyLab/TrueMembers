from __future__ import annotations

from enum import Enum


class DetectionMethod(str, Enum):
    JACKHMMER = "JACKHMMER"
    FOLDSEEK = "FOLDSEEK"
    JACKHMMER_AND_FOLDSEEK = "JACKHMMER_AND_FOLDSEEK"

    def merged(self, other: DetectionMethod) -> DetectionMethod:
        if self is other:
            return self
        return DetectionMethod.JACKHMMER_AND_FOLDSEEK
