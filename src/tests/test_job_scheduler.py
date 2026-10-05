"""Unit tests for JobScheduler.py — convergence and failed-job retry using DB."""

import sqlite3
from unittest.mock import patch

import pytest

from src.scripts.JobStatus import JobStatus
from src.scripts.pipeline_db import init_db, insert_used_seed


def _make_scheduler(db_path):
    from src.JobScheduler import JobScheduler

    class FakeScheduler(JobScheduler):
        def __init__(self, db_path):
            self.db_path = db_path
            self.current_iteration = 1

    return FakeScheduler(db_path)


class TestComputeConvergence:
    def _make_db(self, tmp_path, rows):
        db_path = tmp_path / "test.db"
        init_db(db_path)
        with sqlite3.connect(str(db_path)) as conn:
            for uid, uniprot, iteration, new_members in rows:
                conn.execute(
                    "INSERT INTO used_seeds (unique_id, uniprot_id, iteration, new_members, failed) "
                    "VALUES (?, ?, ?, ?, 0)",
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


class TestRestartFailedJobs:
    def _make_seed(self, unique_id, uniprot_id):
        from pathlib import Path as _Path
        from src.Hit import Hit, DetectionMethod
        from src.Seed import Seed

        hit = Hit(
            unique_id=unique_id,
            uniprot_id=uniprot_id,
            region_start=1,
            region_end=10,
            sequence_path=_Path("/tmp/x.fasta"),
            structure_path=_Path("/tmp/x.cif"),
            detection_method=DetectionMethod.BOTH,
            qtmscore=0.8,
            is_member=True,
        )
        return Seed(hit=hit, job_config_path=_Path("/tmp/job.yaml"),
                    slurm_script_path=_Path("/tmp/submit.sh"))

    def test_first_failure_restarts_job(self, tmp_path):
        db_path = tmp_path / "test.db"
        init_db(db_path)
        insert_used_seed(db_path, "uid1", "A", 1)

        s = _make_scheduler(db_path)
        s.jobs = {10: {"status": JobStatus.FAILED,
                       "seed": self._make_seed("uid1", "A"),
                       "iteration": 1, "retried": False}}
        with patch("subprocess.run") as mock_run:
            mock_run.return_value.stdout = "Submitted batch job 42\n"
            s.restart_failed_jobs()

        assert s.jobs[10]["status"] == JobStatus.HARVESTED
        assert 42 in s.jobs
        assert s.jobs[42]["status"] == JobStatus.PENDING
        assert s.jobs[42]["retried"] is True

    def test_second_failure_drops_job_and_marks_seed(self, tmp_path):
        db_path = tmp_path / "test.db"
        init_db(db_path)
        insert_used_seed(db_path, "uid1", "A", 1)

        s = _make_scheduler(db_path)
        s.jobs = {10: {"status": JobStatus.FAILED,
                       "seed": self._make_seed("uid1", "A"),
                       "iteration": 2, "retried": True}}
        s.restart_failed_jobs()

        assert s.jobs[10]["status"] == JobStatus.DROPPED
        with sqlite3.connect(str(db_path)) as conn:
            row = conn.execute(
                "SELECT failed, new_members FROM used_seeds WHERE unique_id = 'uid1'"
            ).fetchone()
        assert row[0] == 1
        assert row[1] == 0

    def test_no_failed_jobs_does_nothing(self, tmp_path):
        db_path = tmp_path / "test.db"
        init_db(db_path)
        s = _make_scheduler(db_path)
        s.jobs = {}
        s.restart_failed_jobs()

    def test_already_harvested_jobs_ignored(self, tmp_path):
        db_path = tmp_path / "test.db"
        init_db(db_path)
        insert_used_seed(db_path, "uid1", "A", 1)

        s = _make_scheduler(db_path)
        s.jobs = {10: {"status": JobStatus.HARVESTED,
                       "seed": self._make_seed("uid1", "A"),
                       "iteration": 1, "retried": False}}
        s.restart_failed_jobs()

        with sqlite3.connect(str(db_path)) as conn:
            row = conn.execute(
                "SELECT failed FROM used_seeds WHERE unique_id = 'uid1'"
            ).fetchone()
        assert row[0] == 0


class TestGetJobs:
    def test_returns_matching_job_ids(self, tmp_path):
        db_path = tmp_path / "test.db"
        init_db(db_path)
        s = _make_scheduler(db_path)
        s.jobs = {
            1: {"status": JobStatus.RUNNING},
            2: {"status": JobStatus.FAILED},
            3: {"status": JobStatus.FAILED},
            4: {"status": JobStatus.DROPPED},
        }
        assert sorted(s.get_jobs(JobStatus.FAILED)) == [2, 3]
        assert s.get_jobs(JobStatus.DROPPED) == [4]
        assert s.get_jobs(JobStatus.COMPLETED) == []