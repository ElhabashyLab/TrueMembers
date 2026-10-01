"""Unit tests for JobScheduler.py — convergence computation and seed selection."""

import pandas as pd
import pytest

from src.scripts.JobStatus import JobStatus


class TestComputeConvergence:
    @pytest.fixture
    def scheduler(self, tmp_path):
        from src.JobScheduler import JobScheduler

        class FakeScheduler(JobScheduler):
            def __init__(self, used_seeds_path):
                self.used_seeds_path = used_seeds_path

        return FakeScheduler

    def test_converged_when_all_below_floor(self, scheduler, tmp_path):
        used_seeds = tmp_path / "used_seeds.csv"
        df = pd.DataFrame({
            "id": ["A", "B", "C", "D", "E"],
            "region_start": [10, 20, 30, 40, 50],
            "region_end": [100, 100, 100, 100, 100],
            "iteration": [1, 2, 3, 4, 5],
            "new_true_members": [2, 3, 1, 0, 4],
        })
        df.to_csv(used_seeds, index=False)
        s = scheduler(str(used_seeds))
        converged, count = s.compute_convergence(absolute_floor=5, patience_window=5)
        assert converged is True
        assert count == 5

    def test_not_converged_when_one_above_floor(self, scheduler, tmp_path):
        used_seeds = tmp_path / "used_seeds.csv"
        df = pd.DataFrame({
            "id": ["A", "B", "C", "D", "E"],
            "region_start": [10, 20, 30, 40, 50],
            "region_end": [100, 100, 100, 100, 100],
            "iteration": [1, 2, 3, 4, 5],
            "new_true_members": [2, 3, 1, 8, 4],
        })
        df.to_csv(used_seeds, index=False)
        s = scheduler(str(used_seeds))
        converged, count = s.compute_convergence(absolute_floor=5, patience_window=5)
        assert converged is False
        assert count == 4

    def test_too_few_jobs_for_window(self, scheduler, tmp_path):
        used_seeds = tmp_path / "used_seeds.csv"
        df = pd.DataFrame({
            "id": ["A", "B"],
            "region_start": [10, 20],
            "region_end": [100, 100],
            "iteration": [1, 2],
            "new_true_members": [2, 3],
        })
        df.to_csv(used_seeds, index=False)
        s = scheduler(str(used_seeds))
        converged, count = s.compute_convergence(absolute_floor=5, patience_window=5)
        assert converged is False
        assert count == 2

    def test_ignores_unfinished_jobs_with_na(self, scheduler, tmp_path):
        used_seeds = tmp_path / "used_seeds.csv"
        df = pd.DataFrame({
            "id": ["A", "B", "C", "D", "E", "F"],
            "region_start": [10, 20, 30, 40, 50, 60],
            "region_end": [100, 100, 100, 100, 100, 100],
            "iteration": [1, 2, 3, 4, 5, 6],
            "new_true_members": [2, 3, 1, 0, 4, pd.NA],
        })
        df.to_csv(used_seeds, index=False)
        s = scheduler(str(used_seeds))
        converged, count = s.compute_convergence(absolute_floor=5, patience_window=5)
        assert converged is True
        assert count == 5

    def test_exactly_at_floor_not_below(self, scheduler, tmp_path):
        used_seeds = tmp_path / "used_seeds.csv"
        df = pd.DataFrame({
            "id": ["A", "B", "C"],
            "region_start": [10, 20, 30],
            "region_end": [100, 100, 100],
            "iteration": [1, 2, 3],
            "new_true_members": [5, 5, 5],
        })
        df.to_csv(used_seeds, index=False)
        s = scheduler(str(used_seeds))
        converged, _ = s.compute_convergence(absolute_floor=5, patience_window=3)
        assert converged is False

    def test_empty_used_seeds(self, scheduler, tmp_path):
        used_seeds = tmp_path / "used_seeds.csv"
        df = pd.DataFrame(columns=["id", "region_start", "region_end", "iteration", "new_true_members"])
        df.to_csv(used_seeds, index=False)
        s = scheduler(str(used_seeds))
        converged, count = s.compute_convergence(absolute_floor=5, patience_window=5)
        assert converged is False
        assert count == 0


class TestGetHighestScoringSeed:
    def test_picks_highest_qtmscore_unused(self, tmp_path):
        from src.JobScheduler import JobScheduler

        class FakeScheduler(JobScheduler):
            def __init__(self, seed_list_path, used_seeds_path):
                self.members_path = seed_list_path
                self.used_seeds_path = used_seeds_path

        seeds_path = tmp_path / "seeds.csv"
        used_path = tmp_path / "used_seeds.csv"

        pd.DataFrame({
            "id": ["A", "B", "C"],
            "qtmscore": [0.6, 0.9, 0.75],
            "region_start": [10, 20, 30],
            "region_end": [100, 100, 100],
        }).to_csv(seeds_path, index=False)

        pd.DataFrame({
            "id": ["B"],
            "region_start": [20],
            "region_end": [100],
            "iteration": [1],
            "new_true_members": [5],
        }).to_csv(used_path, index=False)

        s = FakeScheduler(str(seeds_path), str(used_path))
        result = s.get_highest_scoring_seed()
        assert result.iloc[0]["id"] == "C"

    def test_picks_only_available_when_top_is_used(self, tmp_path):
        from src.JobScheduler import JobScheduler

        class FakeScheduler(JobScheduler):
            def __init__(self, seed_list_path, used_seeds_path):
                self.members_path = seed_list_path
                self.used_seeds_path = used_seeds_path

        seeds_path = tmp_path / "seeds.csv"
        used_path = tmp_path / "used_seeds.csv"

        pd.DataFrame({
            "id": ["A", "B"],
            "qtmscore": [0.9, 0.6],
            "region_start": [10, 20],
            "region_end": [100, 100],
        }).to_csv(seeds_path, index=False)

        pd.DataFrame({
            "id": ["A"],
            "region_start": [10],
            "region_end": [100],
            "iteration": [1],
            "new_true_members": [3],
        }).to_csv(used_path, index=False)

        s = FakeScheduler(str(seeds_path), str(used_path))
        result = s.get_highest_scoring_seed()
        assert result.iloc[0]["id"] == "B"

    def test_single_seed(self, tmp_path):
        from src.JobScheduler import JobScheduler

        class FakeScheduler(JobScheduler):
            def __init__(self, seed_list_path, used_seeds_path):
                self.members_path = seed_list_path
                self.used_seeds_path = used_seeds_path

        seeds_path = tmp_path / "seeds.csv"
        used_path = tmp_path / "used_seeds.csv"

        pd.DataFrame({
            "id": ["A"],
            "qtmscore": [0.8],
            "region_start": [10],
            "region_end": [100],
        }).to_csv(seeds_path, index=False)

        pd.DataFrame(columns=["id", "region_start", "region_end", "iteration", "new_true_members"]).to_csv(used_path, index=False)

        s = FakeScheduler(str(seeds_path), str(used_path))
        result = s.get_highest_scoring_seed()
        assert result.iloc[0]["id"] == "A"


class TestRemoveSeedsOfFailedJobs:
    @pytest.fixture
    def setup_scheduler(self, tmp_path):
        from src.JobScheduler import JobScheduler

        used_seeds_path = str(tmp_path / "used_seeds.csv")

        class FakeScheduler(JobScheduler):
            def __init__(self):
                self.used_seeds_path = used_seeds_path
                self._seeds_retried = set()
                self.jobs = {}

        return FakeScheduler, used_seeds_path, tmp_path

    def test_first_failure_frees_seed(self, setup_scheduler, mocker):
        FakeScheduler, used_seeds_path, tmp_path = setup_scheduler

        df = pd.DataFrame({
            "id": ["P12345"],
            "region_start": [10],
            "region_end": [100],
            "iteration": [3],
            "new_true_members": [pd.NA],
        })
        df.to_csv(used_seeds_path, index=False)

        mocker.patch(
            "True_Members.JobScheduler.read_yaml_config",
            return_value={
                "global": {"iteration": 3},
                "run": {"seed_uid": "P12345"},
            },
        )

        s = FakeScheduler()
        s.jobs = {101: {"status": JobStatus.FAILED, "config_path": "/fake/config.yaml"}}

        s.remove_seeds_of_failed_jobs_from_used_seeds()

        result = pd.read_csv(used_seeds_path)
        assert len(result) == 0
        assert "P12345" in s._seeds_retried
        assert s.jobs[101]["status"] == JobStatus.HARVESTED

    def test_second_failure_keeps_seed(self, setup_scheduler, mocker):
        FakeScheduler, used_seeds_path, tmp_path = setup_scheduler

        df = pd.DataFrame({
            "id": ["P12345"],
            "region_start": [10],
            "region_end": [100],
            "iteration": [5],
            "new_true_members": [pd.NA],
        })
        df.to_csv(used_seeds_path, index=False)

        mocker.patch(
            "True_Members.JobScheduler.read_yaml_config",
            return_value={
                "global": {"iteration": 5},
                "run": {"seed_uid": "P12345"},
            },
        )

        s = FakeScheduler()
        s._seeds_retried.add("P12345")
        s.jobs = {102: {"status": JobStatus.FAILED, "config_path": "/fake/config.yaml"}}

        s.remove_seeds_of_failed_jobs_from_used_seeds()

        result = pd.read_csv(used_seeds_path)
        assert len(result) == 1
        assert result.iloc[0]["new_true_members"] == 0
        assert s.jobs[102]["status"] == JobStatus.HARVESTED

    def test_two_different_seeds_each_retried_once(self, setup_scheduler, mocker):
        FakeScheduler, used_seeds_path, tmp_path = setup_scheduler

        df = pd.DataFrame({
            "id": ["A", "B"],
            "region_start": [10, 20],
            "region_end": [100, 200],
            "iteration": [1, 2],
            "new_true_members": [pd.NA, pd.NA],
        })
        df.to_csv(used_seeds_path, index=False)

        configs = {
            "/fake/config_a.yaml": {"global": {"iteration": 1}, "run": {"seed_uid": "A"}},
            "/fake/config_b.yaml": {"global": {"iteration": 2}, "run": {"seed_uid": "B"}},
        }
        mocker.patch(
            "True_Members.JobScheduler.read_yaml_config",
            side_effect=lambda path: configs[path],
        )

        s = FakeScheduler()
        s.jobs = {
            201: {"status": JobStatus.FAILED, "config_path": "/fake/config_a.yaml"},
            202: {"status": JobStatus.FAILED, "config_path": "/fake/config_b.yaml"},
        }

        s.remove_seeds_of_failed_jobs_from_used_seeds()

        result = pd.read_csv(used_seeds_path)
        assert len(result) == 0
        assert "A" in s._seeds_retried
        assert "B" in s._seeds_retried
        assert s.jobs[201]["status"] == JobStatus.HARVESTED
        assert s.jobs[202]["status"] == JobStatus.HARVESTED

    def test_seed_a_fails_twice_seed_b_fails_once(self, setup_scheduler, mocker):
        FakeScheduler, used_seeds_path, tmp_path = setup_scheduler

        df = pd.DataFrame({
            "id": ["A", "B"],
            "region_start": [10, 20],
            "region_end": [100, 200],
            "iteration": [1, 2],
            "new_true_members": [pd.NA, pd.NA],
        })
        df.to_csv(used_seeds_path, index=False)

        configs = {
            "/fake/config_a.yaml": {"global": {"iteration": 1}, "run": {"seed_uid": "A"}},
            "/fake/config_b.yaml": {"global": {"iteration": 2}, "run": {"seed_uid": "B"}},
        }
        mocker.patch(
            "True_Members.JobScheduler.read_yaml_config",
            side_effect=lambda path: configs[path],
        )

        s = FakeScheduler()
        s._seeds_retried.add("A")
        s.jobs = {
            301: {"status": JobStatus.FAILED, "config_path": "/fake/config_a.yaml"},
            302: {"status": JobStatus.FAILED, "config_path": "/fake/config_b.yaml"},
        }

        s.remove_seeds_of_failed_jobs_from_used_seeds()

        result = pd.read_csv(used_seeds_path)
        assert len(result) == 1
        assert result.iloc[0]["id"] == "A"
        assert result.iloc[0]["new_true_members"] == 0
        assert s.jobs[301]["status"] == JobStatus.HARVESTED
        assert s.jobs[302]["status"] == JobStatus.HARVESTED

    def test_no_failed_jobs_does_nothing(self, setup_scheduler, mocker):
        FakeScheduler, used_seeds_path, tmp_path = setup_scheduler

        df = pd.DataFrame({
            "id": ["P12345"],
            "region_start": [10],
            "region_end": [100],
            "iteration": [1],
            "new_true_members": [pd.NA],
        })
        df.to_csv(used_seeds_path, index=False)

        mock_read = mocker.patch("True_Members.JobScheduler.read_yaml_config")

        s = FakeScheduler()
        s.jobs = {401: {"status": JobStatus.COMPLETED, "config_path": "/fake/config.yaml"}}

        s.remove_seeds_of_failed_jobs_from_used_seeds()

        result = pd.read_csv(used_seeds_path)
        assert len(result) == 1
        mock_read.assert_not_called()

    def test_config_read_error_does_not_crash(self, setup_scheduler, mocker):
        FakeScheduler, used_seeds_path, tmp_path = setup_scheduler

        df = pd.DataFrame({
            "id": ["P12345"],
            "region_start": [10],
            "region_end": [100],
            "iteration": [3],
            "new_true_members": [pd.NA],
        })
        df.to_csv(used_seeds_path, index=False)

        mocker.patch(
            "True_Members.JobScheduler.read_yaml_config",
            side_effect=OSError("Config file missing"),
        )

        s = FakeScheduler()
        s.jobs = {501: {"status": JobStatus.FAILED, "config_path": "/fake/config.yaml"}}

        s.remove_seeds_of_failed_jobs_from_used_seeds()

        result = pd.read_csv(used_seeds_path)
        assert len(result) == 1
        assert s.jobs[501]["status"] == JobStatus.HARVESTED