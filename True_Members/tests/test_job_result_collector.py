"""Unit tests for scripts/JobResultCollector.py — merge logic for global state."""

import pandas as pd
import pytest

from True_Members.scripts.JobResultCollector import JobResultCollector


class TestMergeLocalTrueMembers:
    def test_new_ids_appended(self, tmp_path):
        seeds_path = tmp_path / "seeds.csv"
        pd.DataFrame({
            "id": ["P1", "P2"],
            "qtmscore": [0.8, 0.7],
            "region_start": [10, 20],
            "region_end": [100, 100],
            "DomainA": [True, True],
        }).to_csv(seeds_path, index=False)

        local_results = tmp_path / "local_results"
        local_results.mkdir()
        pd.DataFrame({
            "id": ["P3"],
            "qtmscore": [0.9],
            "region_start": [30],
            "region_end": [100],
            "DomainA": [True],
        }).to_csv(local_results / "True_Members.csv", index=False)

        collector = JobResultCollector(str(seeds_path), str(tmp_path / "used.csv"), str(tmp_path / "all.csv"))
        merged, new_count = collector.merge_local_true_members_with_global_seeds(1, str(local_results))

        assert len(merged) == 3
        assert new_count == 1
        assert "P3" in merged["id"].values

    def test_duplicate_ids_merged_with_median_qtmscore(self, tmp_path):
        seeds_path = tmp_path / "seeds.csv"
        pd.DataFrame({
            "id": ["P1"],
            "qtmscore": [0.8],
            "region_start": [10],
            "region_end": [100],
            "DomainA": [True],
        }).to_csv(seeds_path, index=False)

        local_results = tmp_path / "local_results"
        local_results.mkdir()
        pd.DataFrame({
            "id": ["P1"],
            "qtmscore": [0.9],
            "region_start": [10],
            "region_end": [100],
            "DomainA": [True],
        }).to_csv(local_results / "True_Members.csv", index=False)

        collector = JobResultCollector(str(seeds_path), str(tmp_path / "used.csv"), str(tmp_path / "all.csv"))
        merged, new_count = collector.merge_local_true_members_with_global_seeds(1, str(local_results))

        assert len(merged) == 1
        assert new_count == 0
        assert merged.iloc[0]["qtmscore"] == pytest.approx(0.85)

    def test_region_expanded_on_merge(self, tmp_path):
        seeds_path = tmp_path / "seeds.csv"
        pd.DataFrame({
            "id": ["P1"],
            "qtmscore": [0.8],
            "region_start": [20],
            "region_end": [80],
            "DomainA": [True],
        }).to_csv(seeds_path, index=False)

        local_results = tmp_path / "local_results"
        local_results.mkdir()
        pd.DataFrame({
            "id": ["P1"],
            "qtmscore": [0.9],
            "region_start": [10],
            "region_end": [100],
            "DomainA": [True],
        }).to_csv(local_results / "True_Members.csv", index=False)

        collector = JobResultCollector(str(seeds_path), str(tmp_path / "used.csv"), str(tmp_path / "all.csv"))
        merged, _ = collector.merge_local_true_members_with_global_seeds(1, str(local_results))

        assert merged.iloc[0]["region_start"] == 10
        assert merged.iloc[0]["region_end"] == 100

    def test_sorted_by_qtmscore_descending(self, tmp_path):
        seeds_path = tmp_path / "seeds.csv"
        pd.DataFrame({
            "id": ["P1"],
            "qtmscore": [0.6],
            "region_start": [10],
            "region_end": [100],
            "DomainA": [True],
        }).to_csv(seeds_path, index=False)

        local_results = tmp_path / "local_results"
        local_results.mkdir()
        pd.DataFrame({
            "id": ["P2"],
            "qtmscore": [0.95],
            "region_start": [20],
            "region_end": [100],
            "DomainA": [True],
        }).to_csv(local_results / "True_Members.csv", index=False)

        collector = JobResultCollector(str(seeds_path), str(tmp_path / "used.csv"), str(tmp_path / "all.csv"))
        merged, _ = collector.merge_local_true_members_with_global_seeds(1, str(local_results))

        assert merged.iloc[0]["id"] == "P2"
        assert merged.iloc[1]["id"] == "P1"


class TestMergeLocalAllMembers:
    def test_any_true_for_hmm_columns(self, tmp_path):
        all_path = tmp_path / "all.csv"
        pd.DataFrame({
            "id": ["P1"],
            "qtmscore": [0.8],
            "region_start": [10],
            "region_end": [100],
            "DomainA": [True],
            "DomainB": [False],
        }).to_csv(all_path, index=False)

        local_results = tmp_path / "local_results"
        local_results.mkdir()
        pd.DataFrame({
            "id": ["P1"],
            "qtmscore": [0.9],
            "region_start": [10],
            "region_end": [100],
            "DomainA": [False],
            "DomainB": [True],
        }).to_csv(local_results / "All_Members.csv", index=False)

        collector = JobResultCollector(str(tmp_path / "seeds.csv"), str(tmp_path / "used.csv"), str(all_path))
        merged = collector.merge_local_all_members_with_global_all_members(1, str(local_results))

        assert merged.iloc[0]["DomainA"] == True
        assert merged.iloc[0]["DomainB"] == True

    def test_appends_new_ids(self, tmp_path):
        all_path = tmp_path / "all.csv"
        pd.DataFrame({
            "id": ["P1"],
            "qtmscore": [0.8],
            "region_start": [10],
            "region_end": [100],
            "DomainA": [True],
        }).to_csv(all_path, index=False)

        local_results = tmp_path / "local_results"
        local_results.mkdir()
        pd.DataFrame({
            "id": ["P2"],
            "qtmscore": [0.7],
            "region_start": [20],
            "region_end": [100],
            "DomainA": [True],
        }).to_csv(local_results / "All_Members.csv", index=False)

        collector = JobResultCollector(str(tmp_path / "seeds.csv"), str(tmp_path / "used.csv"), str(all_path))
        merged = collector.merge_local_all_members_with_global_all_members(1, str(local_results))

        assert len(merged) == 2
        assert set(merged["id"]) == {"P1", "P2"}


class TestAddCurrentSeedInfo:
    def test_updates_new_true_members(self, tmp_path):
        used_path = tmp_path / "used.csv"
        pd.DataFrame({
            "id": ["A", "B"],
            "region_start": [10, 20],
            "region_end": [100, 100],
            "iteration": [1, 2],
            "new_true_members": [pd.NA, pd.NA],
        }).to_csv(used_path, index=False)

        collector = JobResultCollector(str(tmp_path / "seeds.csv"), str(used_path), str(tmp_path / "all.csv"))
        result = collector.add_current_seed_info_to_used_seeds("A", 42)

        assert result.loc[result["id"] == "A", "new_true_members"].iloc[0] == 42
        assert pd.isna(result.loc[result["id"] == "B", "new_true_members"].iloc[0])