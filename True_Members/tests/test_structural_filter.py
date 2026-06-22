"""Unit tests for scripts/structural_filter.py — output parsing and filtering."""

import pandas as pd

from True_Members.scripts.structural_filter import filter_foldseek_output, filter_usalign_output


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