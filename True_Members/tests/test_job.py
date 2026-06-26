"""Unit tests for Job.py — combine_results logic."""

import pandas as pd
import pytest

from True_Members.scripts.types import JobStats


class FakeJob:
    def __init__(self):
        self.stats = JobStats()

    def combine_results(self, uids_with_region, seed_containing_sequences, hmm_search_results, hmms_with_ga_cutoffs):
        combined_results = uids_with_region.merge(seed_containing_sequences, on="id", how="left")
        combined_results = combined_results.dropna(subset=["qtmscore"])
        combined_results = pd.merge(combined_results, hmm_search_results, on="id", how="left")
        expected_hmms = list(hmms_with_ga_cutoffs.keys())
        for col in expected_hmms:
            if col not in combined_results.columns:
                combined_results[col] = False
        combined_results[expected_hmms] = combined_results[expected_hmms].fillna(False).astype(bool)

        all_members = combined_results
        self.stats.all_members = len(all_members)
        true_members = all_members[all_members[expected_hmms].all(axis=1)]
        self.stats.true_members = len(true_members)
        return all_members, true_members


class TestCombineResults:
    def test_all_pass_both_filters(self):
        uids = pd.DataFrame({"id": ["P1", "P2"], "region_start": [10, 20], "region_end": [30, 40]})
        structural = pd.DataFrame({"id": ["P1", "P2"], "qtmscore": [0.85, 0.75]})
        hmm_results = pd.DataFrame({
            "id": ["P1", "P2"],
            "DomainA": [True, True],
            "DomainB": [True, True],
        })
        hmm_cutoffs = {"DomainA": 25.0, "DomainB": 25.0}

        job = FakeJob()
        all_members, true_members = job.combine_results(uids, structural, hmm_results, hmm_cutoffs)

        assert job.stats.all_members == 2
        assert job.stats.true_members == 2
        assert len(true_members) == 2

    def test_structural_filter_drops_sequences(self):
        uids = pd.DataFrame({"id": ["P1", "P2", "P3"], "region_start": [10, 20, 30], "region_end": [30, 40, 50]})
        structural = pd.DataFrame({"id": ["P1", "P3"], "qtmscore": [0.85, 0.75]})
        hmm_results = pd.DataFrame({
            "id": ["P1", "P3"],
            "DomainA": [True, True],
            "DomainB": [True, True],
        })
        hmm_cutoffs = {"DomainA": 25.0, "DomainB": 25.0}

        job = FakeJob()
        all_members, true_members = job.combine_results(uids, structural, hmm_results, hmm_cutoffs)

        assert job.stats.all_members == 2
        assert job.stats.true_members == 2
        assert "P2" not in all_members["id"].values

    def test_missing_domain_columns_filled_with_false(self):
        uids = pd.DataFrame({"id": ["P1", "P2"], "region_start": [10, 20], "region_end": [30, 40]})
        structural = pd.DataFrame({"id": ["P1", "P2"], "qtmscore": [0.85, 0.75]})
        hmm_results = pd.DataFrame({
            "id": ["P1", "P2"],
            "DomainA": [True, False],
        })
        hmm_cutoffs = {"DomainA": 25.0, "DomainB": 25.0}

        job = FakeJob()
        all_members, true_members = job.combine_results(uids, structural, hmm_results, hmm_cutoffs)

        assert "DomainB" in all_members.columns
        assert all_members["DomainB"].tolist() == [False, False]
        assert len(true_members) == 0

    def test_only_sequences_with_all_domains_are_true_members(self):
        uids = pd.DataFrame({"id": ["P1", "P2", "P3"], "region_start": [10, 20, 30], "region_end": [30, 40, 50]})
        structural = pd.DataFrame({"id": ["P1", "P2", "P3"], "qtmscore": [0.85, 0.75, 0.90]})
        hmm_results = pd.DataFrame({
            "id": ["P1", "P2", "P3"],
            "DomainA": [True, True, False],
            "DomainB": [True, False, True],
        })
        hmm_cutoffs = {"DomainA": 25.0, "DomainB": 25.0}

        job = FakeJob()
        all_members, true_members = job.combine_results(uids, structural, hmm_results, hmm_cutoffs)

        assert job.stats.all_members == 3
        assert job.stats.true_members == 1
        assert true_members.iloc[0]["id"] == "P1"

    def test_empty_hmm_results(self):
        uids = pd.DataFrame({"id": ["P1"], "region_start": [10], "region_end": [30]})
        structural = pd.DataFrame({"id": ["P1"], "qtmscore": [0.85]})
        hmm_results = pd.DataFrame(columns=["id"])
        hmm_cutoffs = {"DomainA": 25.0}

        job = FakeJob()
        all_members, true_members = job.combine_results(uids, structural, hmm_results, hmm_cutoffs)

        assert job.stats.all_members == 1
        assert "DomainA" in all_members.columns
        assert job.stats.true_members == 0

    def test_empty_uids(self):
        uids = pd.DataFrame(columns=["id", "region_start", "region_end"])
        structural = pd.DataFrame(columns=["id", "qtmscore"])
        hmm_results = pd.DataFrame(columns=["id"])
        hmm_cutoffs = {"DomainA": 25.0}

        job = FakeJob()
        all_members, true_members = job.combine_results(uids, structural, hmm_results, hmm_cutoffs)

        assert job.stats.all_members == 0
        assert job.stats.true_members == 0
        assert all_members.empty

    def test_empty_structural_filter_with_uids(self):
        uids = pd.DataFrame({"id": ["P1", "P2"], "region_start": [10, 20], "region_end": [30, 40]})
        structural = pd.DataFrame(columns=["id", "qtmscore"])
        hmm_results = pd.DataFrame(columns=["id"])
        hmm_cutoffs = {"DomainA": 25.0}

        job = FakeJob()
        all_members, true_members = job.combine_results(uids, structural, hmm_results, hmm_cutoffs)

        assert job.stats.all_members == 0
        assert job.stats.true_members == 0
        assert all_members.empty

    def test_empty_hmm_results_returns_no_true_members(self):
        uids = pd.DataFrame({"id": ["P1"], "region_start": [10], "region_end": [30]})
        structural = pd.DataFrame({"id": ["P1"], "qtmscore": [0.85]})
        hmm_results = pd.DataFrame(columns=["id"])
        hmm_cutoffs = {"DomainA": 25.0}

        job = FakeJob()
        all_members, true_members = job.combine_results(uids, structural, hmm_results, hmm_cutoffs)

        assert job.stats.all_members == 1
        assert "DomainA" in all_members.columns
        assert all_members["DomainA"].tolist() == [False]
        assert job.stats.true_members == 0