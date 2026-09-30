"""Unit tests for scripts/get_structures_and_sequences.py — annotation file handling and ID merging."""

import os

import pandas as pd

from src.scripts.get_structures_and_sequences import (
    combine_foldseek_and_jackhmmer_ids,
    get_annotation_files,
    get_jackhmmer_and_foldseek_ids_with_regions,
    get_unique_ids_with_regions_from_annotation_files,
)
from src.scripts.types import JobStats


class TestGetAnnotationFiles:
    def test_finds_annotation_files(self, tmp_path):
        seed = ["Q5SK37", 507, 569]
        bitscore = 0.8
        prefix = f"Q5SK37_507-569_b{bitscore}"
        align_dir = os.path.join(str(tmp_path), "output", prefix, "align")
        os.makedirs(align_dir, exist_ok=True)
        open(os.path.join(align_dir, f"{prefix}_annotation.csv"), "w").write("id\nQ5SK37/507-569")

        result = get_annotation_files(str(tmp_path), seed, [bitscore])
        assert len(result) == 1
        assert "annotation.csv" in result[0]

    def test_missing_annotation_file_warns(self, tmp_path):
        seed = ["Q5SK37", 507, 569]
        result = get_annotation_files(str(tmp_path), seed, [0.2])
        assert result == []

    def test_multiple_bitscores(self, tmp_path):
        seed = ["Q5SK37", 507, 569]
        for b in [0.6, 0.8]:
            prefix = f"Q5SK37_507-569_b{b}"
            align_dir = os.path.join(str(tmp_path), "output", prefix, "align")
            os.makedirs(align_dir, exist_ok=True)
            open(os.path.join(align_dir, f"{prefix}_annotation.csv"), "w").write("id\n")

        result = get_annotation_files(str(tmp_path), seed, [0.6, 0.8])
        assert len(result) == 2


class TestGetUniqueIdsWithRegions:
    def test_extracts_ids_from_real_annotation(self, annotation_csv_small):
        df = get_unique_ids_with_regions_from_annotation_files([annotation_csv_small])
        assert "id" in df.columns
        assert "region_start" in df.columns
        assert "region_end" in df.columns
        assert len(df) > 0
        assert pd.api.types.is_string_dtype(df["id"])

    def test_handles_empty_file_list(self):
        df = get_unique_ids_with_regions_from_annotation_files([])
        assert df.empty
        assert list(df.columns) == ["id", "region_start", "region_end"]

    def test_deduplicates_same_id_across_files(self, tmp_path):
        csv1 = tmp_path / "file1.csv"
        csv2 = tmp_path / "file2.csv"
        csv1_data = pd.DataFrame({"id": ["tr|P12345|PROT/10-20"]})
        csv2_data = pd.DataFrame({"id": ["tr|P12345|PROT/15-25"]})
        csv1_data.to_csv(csv1, index=False)
        csv2_data.to_csv(csv2, index=False)

        df = get_unique_ids_with_regions_from_annotation_files([str(csv1), str(csv2)])
        assert len(df) == 1
        assert df.iloc[0]["id"] == "P12345"
        assert df.iloc[0]["region_start"] == 10
        assert df.iloc[0]["region_end"] == 25

    def test_merges_different_ids_from_multiple_files(self, tmp_path):
        csv1 = tmp_path / "file1.csv"
        csv2 = tmp_path / "file2.csv"
        csv1_data = pd.DataFrame({"id": ["tr|P11111|PROT/10-20"]})
        csv2_data = pd.DataFrame({"id": ["tr|P22222|PROT/30-40"]})
        csv1_data.to_csv(csv1, index=False)
        csv2_data.to_csv(csv2, index=False)

        df = get_unique_ids_with_regions_from_annotation_files([str(csv1), str(csv2)])
        assert len(df) == 2


class TestCombineFoldseekAndJackhmmerIds:
    def test_merges_non_overlapping_ids(self):
        jackhmmer = pd.DataFrame({"id": ["P12345", "Q67890"], "region_start": [10, 30], "region_end": [20, 40]})
        foldseek = pd.DataFrame({"id": ["A11111"], "region_start": [5], "region_end": [15]})

        result = combine_foldseek_and_jackhmmer_ids(jackhmmer, foldseek)
        assert len(result) == 3
        assert set(result["id"]) == {"P12345", "Q67890", "A11111"}

    def test_merges_overlapping_ids_with_averaged_regions(self):
        jackhmmer = pd.DataFrame({"id": ["P12345"], "region_start": [10], "region_end": [20]})
        foldseek = pd.DataFrame({"id": ["P12345"], "region_start": [14], "region_end": [24]})

        result = combine_foldseek_and_jackhmmer_ids(jackhmmer, foldseek)
        assert len(result) == 1
        assert result.iloc[0]["region_start"] == 12
        assert result.iloc[0]["region_end"] == 22

    def test_empty_foldseek_df(self):
        jackhmmer = pd.DataFrame({"id": ["P12345"], "region_start": [10], "region_end": [20]})
        foldseek = pd.DataFrame({"id": pd.Series([], dtype=str), "region_start": pd.Series([], dtype="Int64"), "region_end": pd.Series([], dtype="Int64")})

        result = combine_foldseek_and_jackhmmer_ids(jackhmmer, foldseek)
        assert len(result) == 1
        assert result.iloc[0]["id"] == "P12345"

    def test_both_empty(self):
        jackhmmer = pd.DataFrame({"id": pd.Series([], dtype=str), "region_start": pd.Series([], dtype="Int64"), "region_end": pd.Series([], dtype="Int64")})
        foldseek = jackhmmer.copy()

        result = combine_foldseek_and_jackhmmer_ids(jackhmmer, foldseek)
        assert result.empty


class TestGetJackhmmerAndFoldseekIds:
    def test_missing_foldseek_csv_returns_empty_df(self, tmp_path):
        evcouplings_dir = str(tmp_path / "evcouplings")
        foldseek_dir = str(tmp_path / "foldseek")
        os.makedirs(foldseek_dir, exist_ok=True)

        stats = JobStats()
        seed = ["Q5SK37", 507, 569]
        bit_scores = [0.8]

        jackhmmer_ids, foldseek_ids = get_jackhmmer_and_foldseek_ids_with_regions(
            stats, seed, bit_scores, evcouplings_dir, foldseek_dir
        )

        assert isinstance(jackhmmer_ids, pd.DataFrame)
        assert isinstance(foldseek_ids, pd.DataFrame)
        assert list(foldseek_ids.columns) == ["id", "region_start", "region_end"]
        assert foldseek_ids.empty

    def test_foldseek_csv_exists_but_empty(self, tmp_path):
        evcouplings_dir = str(tmp_path / "evcouplings")
        foldseek_dir = str(tmp_path / "foldseek")
        os.makedirs(foldseek_dir, exist_ok=True)

        empty_csv = pd.DataFrame(columns=["id", "region_start", "region_end"])
        empty_csv.to_csv(f"{foldseek_dir}/foldseek_webserver_output.csv", index=False)

        stats = JobStats()
        seed = ["Q5SK37", 507, 569]
        bit_scores = [0.8]

        jackhmmer_ids, foldseek_ids = get_jackhmmer_and_foldseek_ids_with_regions(
            stats, seed, bit_scores, evcouplings_dir, foldseek_dir
        )

        assert foldseek_ids.empty
        assert list(foldseek_ids.columns) == ["id", "region_start", "region_end"]

    def test_no_annotation_files_all_empty(self, tmp_path):
        evcouplings_dir = str(tmp_path / "evcouplings")
        foldseek_dir = str(tmp_path / "foldseek")
        os.makedirs(foldseek_dir, exist_ok=True)

        stats = JobStats()
        seed = ["Q5SK37", 507, 569]
        bit_scores = [0.2]

        jackhmmer_ids, foldseek_ids = get_jackhmmer_and_foldseek_ids_with_regions(
            stats, seed, bit_scores, evcouplings_dir, foldseek_dir
        )

        assert jackhmmer_ids.empty
        assert stats.annotation_files_jackhmmer == 0
        assert stats.ids_jackhmmer == 0