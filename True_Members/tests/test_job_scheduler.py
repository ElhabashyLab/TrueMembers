"""Unit tests for JobScheduler.py — convergence computation and seed selection."""

import pandas as pd
import pytest

from True_Members.scripts.JobStatus import JobStatus


class TestComputeConvergence:
    @pytest.fixture
    def scheduler(self, tmp_path):
        from True_Members.JobScheduler import JobScheduler

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
        from True_Members.JobScheduler import JobScheduler

        class FakeScheduler(JobScheduler):
            def __init__(self, seed_list_path, used_seeds_path):
                self.seed_list_path = seed_list_path
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
        from True_Members.JobScheduler import JobScheduler

        class FakeScheduler(JobScheduler):
            def __init__(self, seed_list_path, used_seeds_path):
                self.seed_list_path = seed_list_path
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
        from True_Members.JobScheduler import JobScheduler

        class FakeScheduler(JobScheduler):
            def __init__(self, seed_list_path, used_seeds_path):
                self.seed_list_path = seed_list_path
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