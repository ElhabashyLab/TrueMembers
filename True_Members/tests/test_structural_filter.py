"""Unit tests for scripts/structural_filter.py — output parsing and filtering."""

import os

import pandas as pd

from True_Members.scripts.structural_filter import filter_foldseek_output, filter_usalign_output, run_structural_similarity_filter
from True_Members.scripts.types import JobStats


class TestFilterFoldseekOutput:
    def test_filters_by_tm_score(self, tmp_path):
        content = (
            "query1\ttarget1\t0.0\t0.85\t0.80\n"
            "query1\ttarget2\t0.0\t0.55\t0.50\n"
            "query1\ttarget3\t0.0\t0.72\t0.70\n"
        )
        input_file = tmp_path / "foldseek.tsv"
        input_file.write_text(content)

        result = filter_foldseek_output(str(input_file), TM_SCORE_CUTOFF=0.6)
        assert len(result) == 2
        assert set(result["target"]) == {"target1", "target3"}

    def test_no_structures_pass(self, tmp_path):
        content = "query1\ttarget1\t0.0\t0.30\t0.30\n"
        input_file = tmp_path / "foldseek.tsv"
        input_file.write_text(content)

        result = filter_foldseek_output(str(input_file), TM_SCORE_CUTOFF=0.6)
        assert len(result) == 0

    def test_all_structures_pass(self, tmp_path):
        content = "query1\ttarget1\t0.0\t0.95\t0.90\n"
        input_file = tmp_path / "foldseek.tsv"
        input_file.write_text(content)

        result = filter_foldseek_output(str(input_file), TM_SCORE_CUTOFF=0.6)
        assert len(result) == 1

    def test_exact_boundary(self, tmp_path):
        content = "query1\ttarget1\t0.0\t0.6\t0.6\n"
        input_file = tmp_path / "foldseek.tsv"
        input_file.write_text(content)

        result = filter_foldseek_output(str(input_file), TM_SCORE_CUTOFF=0.6)
        assert len(result) == 1


class TestFilterUsalignOutput:
    def test_filters_by_tm_score(self, tmp_path):
        content = (
            "#PDBchain1\tPDBchain2\tTM1\tTM2\tRMSD\tID1\tID2\tIDfull\tL1\tL2\tLali\n"
            "target1\tquery1\t0.50\t0.88\t2.5\t0.45\t0.45\t0.45\t100\t100\t95\n"
            "target2\tquery1\t0.45\t0.52\t3.1\t0.30\t0.30\t0.30\t100\t100\t90\n"
            "target3\tquery1\t0.60\t0.75\t1.8\t0.55\t0.55\t0.55\t100\t100\t97\n"
        )
        input_file = tmp_path / "usalign.txt"
        input_file.write_text(content)

        result = filter_usalign_output(str(input_file), TM_SCORE_CUTOFF=0.6)
        assert len(result) == 2
        assert set(result["target"]) == {"target1", "target3"}

    def test_skips_slurm_footer(self, tmp_path):
        content = (
            "#PDBchain1\tPDBchain2\tTM1\tTM2\tRMSD\tID1\tID2\tIDfull\tL1\tL2\tLali\n"
            "target1\tquery1\t0.50\t0.85\t2.5\t0.45\t0.45\t0.45\t100\t100\t95\n"
            "== Job statistics ==\n"
            "== Some slurm output ==\n"
        )
        input_file = tmp_path / "usalign.txt"
        input_file.write_text(content)

        result = filter_usalign_output(str(input_file), TM_SCORE_CUTOFF=0.6)
        assert len(result) == 1

    def test_no_structures_pass(self, tmp_path):
        content = (
            "#PDBchain1\tPDBchain2\tTM1\tTM2\tRMSD\tID1\tID2\tIDfull\tL1\tL2\tLali\n"
            "target1\tquery1\t0.10\t0.15\t5.0\t0.05\t0.05\t0.05\t100\t100\t50\n"
        )
        input_file = tmp_path / "usalign.txt"
        input_file.write_text(content)

        result = filter_usalign_output(str(input_file), TM_SCORE_CUTOFF=0.6)
        assert len(result) == 0

    def test_cleans_path_extensions_from_ids(self, tmp_path):
        content = (
            "#PDBchain1\tPDBchain2\tTM1\tTM2\tRMSD\tID1\tID2\tIDfull\tL1\tL2\tLali\n"
            "target_dir/P12345.cif\tquery_dir/query.cif\t0.50\t0.85\t2.5\t0.45\t0.45\t0.45\t100\t100\t95\n"
        )
        input_file = tmp_path / "usalign.txt"
        input_file.write_text(content)

        result = filter_usalign_output(str(input_file), TM_SCORE_CUTOFF=0.6)
        assert result.iloc[0]["target"] == "P12345"
        assert result.iloc[0]["query"] == "query"


class TestRunStructuralSimilarityFilter:
    def test_empty_uids_with_region_returns_empty_df(self, mocker, tmp_path):
        structures_dir = str(tmp_path / "structures")
        os.makedirs(structures_dir, exist_ok=True)
        (tmp_path / "structures" / "P12345.cif").write_text("")

        foldseek_dir = str(tmp_path / "foldseek")
        usalign_dir = str(tmp_path / "usalign")
        os.makedirs(foldseek_dir, exist_ok=True)
        os.makedirs(usalign_dir, exist_ok=True)

        foldseek_output = str(tmp_path / "foldseek" / "output_foldseek.csv")
        (tmp_path / "foldseek" / "output_foldseek.csv").write_text(
            "query1\tP12345\t0.0\t0.85\t0.80\n"
        )

        usalign_output = str(tmp_path / "usalign" / "usalign_results.txt")
        (tmp_path / "usalign" / "usalign_results.txt").write_text(
            "#PDBchain1\tPDBchain2\tTM1\tTM2\tRMSD\tID1\tID2\tIDfull\tL1\tL2\tLali\n"
            "P12345\tquery1\t0.50\t0.88\t2.5\t0.45\t0.45\t0.45\t100\t100\t95\n"
        )

        mock_run_foldseek = mocker.patch(
            "True_Members.scripts.structural_filter.run_foldseek",
            return_value=foldseek_output,
        )
        mock_run_usalign = mocker.patch(
            "True_Members.scripts.structural_filter.run_usalign",
            return_value=usalign_output,
        )

        stats = JobStats()
        empty_uids = pd.DataFrame(columns=["id", "region_start", "region_end"])

        result = run_structural_similarity_filter(
            query_structure_roi="/fake/query.cif",
            stats=stats,
            TM_SCORE_CUTOFF=0.6,
            structures_dir=structures_dir,
            foldseek_dir=foldseek_dir,
            usalign_dir=usalign_dir,
            uids_with_region=empty_uids,
            structural_filter_tool="both",
        )

        assert isinstance(result, pd.DataFrame)
        assert result.empty
        assert stats.filtered_structures_total == 0
        mock_run_foldseek.assert_called_once()
        mock_run_usalign.assert_called_once()

    def test_empty_uids_usalign_only(self, mocker, tmp_path):
        structures_dir = str(tmp_path / "structures")
        os.makedirs(structures_dir, exist_ok=True)
        (tmp_path / "structures" / "P12345.cif").write_text("")

        usalign_dir = str(tmp_path / "usalign")
        os.makedirs(usalign_dir, exist_ok=True)

        usalign_output = str(tmp_path / "usalign" / "usalign_results.txt")
        (tmp_path / "usalign" / "usalign_results.txt").write_text(
            "#PDBchain1\tPDBchain2\tTM1\tTM2\tRMSD\tID1\tID2\tIDfull\tL1\tL2\tLali\n"
            "P12345\tquery1\t0.50\t0.88\t2.5\t0.45\t0.45\t0.45\t100\t100\t95\n"
        )

        mock_run_usalign = mocker.patch(
            "True_Members.scripts.structural_filter.run_usalign",
            return_value=usalign_output,
        )

        stats = JobStats()
        empty_uids = pd.DataFrame(columns=["id", "region_start", "region_end"])

        result = run_structural_similarity_filter(
            query_structure_roi="/fake/query.cif",
            stats=stats,
            TM_SCORE_CUTOFF=0.6,
            structures_dir=structures_dir,
            foldseek_dir=str(tmp_path / "foldseek"),
            usalign_dir=usalign_dir,
            uids_with_region=empty_uids,
            structural_filter_tool="usalign",
        )

        assert isinstance(result, pd.DataFrame)
        assert result.empty
        assert stats.filtered_structures_total == 0
        mock_run_usalign.assert_called_once()