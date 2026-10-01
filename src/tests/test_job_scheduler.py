"""Unit tests for JobScheduler.py — convergence and failed-job retry using DB."""

import sqlite3

import pytest

from src.scripts.JobStatus import JobStatus
from src.scripts.pipeline_db import init_db, insert_used_seed


def _make_scheduler(db_path):
    from src.JobScheduler import JobScheduler

    class FakeScheduler(JobScheduler):
        def __init__(self, db_path):
            self.db_path = db_path
            self._seeds_retried = set()

    return FakeScheduler(db_path)


class TestComputeConvergence:
    def _make_db(self, tmp_path, rows):
        db_path = tmp_path / "test.db"
        init_db(db_path)
        with sqlite3.connect(str(db_path)) as conn:
            for uid, uniprot, iteration, new_members in rows:
                conn.execute(
                    "INSERT INTO used_seeds (unique_id, uniprot_id, iteration, new_members) "
                    "VALUES (?, ?, ?, ?)",
                    (uid, uniprot, iteration, new_members),
                )
            conn.commit()
        return db_path

    def test_converged_when_all_below_floor(self, tmp_path):
        db = self._make_db(tmp_path, [
            ("a", "A", 1, 2), ("b", "B", 2, 3), ("c", "C", 3, 1),
            ("d", "D", 4, 0), ("e", "E", 5, 4),
        ])
        s = _make_scheduler(db)
        converged, count = s.compute_convergence(absolute_floor=5, patience_window=5)
        assert converged is True
        assert count == 5

    def test_not_converged_when_one_above_floor(self, tmp_path):
        db = self._make_db(tmp_path, [
            ("a", "A", 1, 2), ("b", "B", 2, 3), ("c", "C", 3, 1),
            ("d", "D", 4, 8), ("e", "E", 5, 4),
        ])
        s = _make_scheduler(db)
        converged, count = s.compute_convergence(absolute_floor=5, patience_window=5)
        assert converged is False
        assert count == 4

    def test_too_few_jobs_for_window(self, tmp_path):
        db = self._make_db(tmp_path, [
            ("a", "A", 1, 2), ("b", "B", 2, 3),
        ])
        s = _make_scheduler(db)
        converged, count = s.compute_convergence(absolute_floor=5, patience_window=5)
        assert converged is False
        assert count == 2

    def test_ignores_unfinished_jobs_with_null(self, tmp_path):
        db = self._make_db(tmp_path, [
            ("a", "A", 1, 2), ("b", "B", 2, 3), ("c", "C", 3, 1),
            ("d", "D", 4, 0), ("e", "E", 5, 4), ("f", "F", 6, None),
        ])
        s = _make_scheduler(db)
        converged, count = s.compute_convergence(absolute_floor=5, patience_window=5)
        assert converged is True
        assert count == 5

    def test_exactly_at_floor_not_below(self, tmp_path):
        db = self._make_db(tmp_path, [
            ("a", "A", 1, 5), ("b", "B", 2, 5), ("c", "C", 3, 5),
        ])
        s = _make_scheduler(db)
        converged, _ = s.compute_convergence(absolute_floor=5, patience_window=3)
        assert converged is False

    def test_empty_used_seeds(self, tmp_path):
        db = self._make_db(tmp_path, [])
        s = _make_scheduler(db)
        converged, count = s.compute_convergence(absolute_floor=5, patience_window=5)
        assert converged is False
        assert count == 0


class TestRemoveSeedsOfFailedJobs:
    def test_first_failure_deletes_row(self, tmp_path):
        import yaml
        db_path = tmp_path / "test.db"
        init_db(db_path)
        insert_used_seed(db_path, "uid1", "A", 1)

        config_path = tmp_path / "config.yaml"
        yaml.safe_dump({"seed": {"unique_id": "uid1"}, "pipeline": {"iteration": 1}},
                       open(config_path, "w"), sort_keys=False)

        s = _make_scheduler(db_path)
        s.jobs = {10: {"status": JobStatus.FAILED, "config_path": config_path}}
        s.remove_seeds_of_failed_jobs_from_used_seeds()

        with sqlite3.connect(str(db_path)) as conn:
            rows = conn.execute("SELECT * FROM used_seeds WHERE unique_id = 'uid1'").fetchall()
        assert len(rows) == 0

    def test_second_failure_marks_new_members_zero(self, tmp_path):
        import yaml
        db_path = tmp_path / "test.db"
        init_db(db_path)
        insert_used_seed(db_path, "uid1", "A", 1)

        config_path = tmp_path / "config.yaml"
        yaml.safe_dump({"seed": {"unique_id": "uid1"}, "pipeline": {"iteration": 1}},
                       open(config_path, "w"), sort_keys=False)

        s = _make_scheduler(db_path)
        s._seeds_retried.add("uid1")
        s.jobs = {10: {"status": JobStatus.FAILED, "config_path": config_path}}
        s.remove_seeds_of_failed_jobs_from_used_seeds()

        with sqlite3.connect(str(db_path)) as conn:
            row = conn.execute(
                "SELECT new_members FROM used_seeds WHERE unique_id = 'uid1' AND iteration = 1"
            ).fetchone()
        assert row[0] == 0

    def test_no_failed_jobs_does_nothing(self, tmp_path):
        db_path = tmp_path / "test.db"
        init_db(db_path)
        s = _make_scheduler(db_path)
        s.jobs = {}
        s.remove_seeds_of_failed_jobs_from_used_seeds()

    def test_config_read_error_does_not_crash(self, tmp_path):
        db_path = tmp_path / "test.db"
        init_db(db_path)

        s = _make_scheduler(db_path)
        s.jobs = {10: {"status": JobStatus.FAILED, "config_path": tmp_path / "nonexistent.yaml"}}
        s.remove_seeds_of_failed_jobs_from_used_seeds()