"""Tests for create_run_summary.py — pipeline run summary generation."""

import json
import os

from create_run_summary import build_summary, find_job_directory, find_used_seeds_file, read_number_tracker


class TestFindUsedSeedsFile:
    def test_finds_csv(self, tmp_path):
        results_dir = tmp_path / "final_true_members_results"
        results_dir.mkdir()
        csv_path = results_dir / "Used_Seeds_i_10.csv"
        csv_path.write_text("id,iteration\n")
        assert find_used_seeds_file(str(tmp_path)) == str(csv_path)

    def test_missing_dir_raises(self, tmp_path):
        import pytest
        with pytest.raises(FileNotFoundError, match="Missing directory"):
            find_used_seeds_file(str(tmp_path))


class TestFindJobDirectory:
    def test_finds_matching_dir(self, tmp_path):
        runs_dir = tmp_path / "runs"
        runs_dir.mkdir()
        job_dir = runs_dir / "job_Q5SK37_iteration_1_2026-06-27_13-33-13"
        job_dir.mkdir(parents=True)
        result = find_job_directory(str(tmp_path), "Q5SK37", 1)
        assert result == str(job_dir)

    def test_no_match_raises(self, tmp_path):
        import pytest
        runs_dir = tmp_path / "runs"
        runs_dir.mkdir()
        with pytest.raises(FileNotFoundError, match="No job directory matching"):
            find_job_directory(str(tmp_path), "NONEXIST", 1)


class TestReadNumberTracker:
    def test_reads_all_stages(self, tmp_path):
        results_dir = tmp_path / "results"
        results_dir.mkdir()
        csv_path = results_dir / "number_tracker.csv"
        csv_path.write_text(
            "Stage,Numbers\n"
            "annotation_files_jackhmmer,4\n"
            "ids_jackhmmer,89882\n"
            "true_members,201\n"
        )
        data = read_number_tracker(str(tmp_path))
        assert data["annotation_files_jackhmmer"] == 4
        assert data["ids_jackhmmer"] == 89882
        assert data["true_members"] == 201

    def test_missing_file_raises(self, tmp_path):
        import pytest
        with pytest.raises(FileNotFoundError, match="Missing number_tracker.csv"):
            read_number_tracker(str(tmp_path))


class TestBuildSummary:
    def test_builds_correct_summary(self, tmp_path):
        results_dir = tmp_path / "final_true_members_results"
        results_dir.mkdir(parents=True)
        used_seeds = results_dir / "Used_Seeds_i_10.csv"
        used_seeds.write_text(
            "id,region_start,region_end,iteration,new_true_members\n"
            "Q5SK37,507,569,1,201.0\n"
        )

        runs_dir = tmp_path / "runs"
        runs_dir.mkdir()
        job_dir = runs_dir / "job_Q5SK37_iteration_1_2026-06-27_13-33-13"
        job_results = job_dir / "results"
        job_results.mkdir(parents=True)
        tracker = job_results / "number_tracker.csv"
        tracker.write_text(
            "Stage,Numbers\n"
            "annotation_files_jackhmmer,4\n"
            "ids_jackhmmer,89882\n"
            "ids_foldseek_webserver,0\n"
            "ids_combined,89882\n"
            "structure_downloads,32325\n"
            "uniprot_records,32325\n"
            "filtered_structures_foldseek,0\n"
            "filtered_structures_usalign,5574\n"
            "filtered_structures_total,5573\n"
            "hmm_search_results,5531\n"
            "all_members,5573\n"
            "true_members,201\n"
        )

        summary = build_summary(str(tmp_path))

        assert "Q5SK37" in summary
        entry = summary["Q5SK37"]
        assert entry["iteration"] == 1
        assert entry["new_true_members"] == 201
        assert entry["annotation_files_jackhmmer"] == 4
        assert entry["ids_jackhmmer"] == 89882
        assert entry["true_members"] == 201

    def test_float_new_true_members(self, tmp_path):
        results_dir = tmp_path / "final_true_members_results"
        results_dir.mkdir(parents=True)
        used_seeds = results_dir / "Used_Seeds_i_10.csv"
        used_seeds.write_text(
            "id,region_start,region_end,iteration,new_true_members\n"
            "Q5SK37,507,569,1,201.0\n"
        )
        runs_dir = tmp_path / "runs"
        runs_dir.mkdir()
        job_dir = runs_dir / "job_Q5SK37_iteration_1_2026-06-27_13-33-13"
        job_results = job_dir / "results"
        job_results.mkdir(parents=True)
        tracker = job_results / "number_tracker.csv"
        tracker.write_text("Stage,Numbers\ntrue_members,201\n")

        summary = build_summary(str(tmp_path))
        assert summary["Q5SK37"]["new_true_members"] == 201