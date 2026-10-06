from datetime import datetime


def duration_seconds(started_at: str | None, finished_at: str | None) -> float | None:
    if not started_at or not finished_at:
        return None
    return (datetime.fromisoformat(finished_at) - datetime.fromisoformat(started_at)).total_seconds()
