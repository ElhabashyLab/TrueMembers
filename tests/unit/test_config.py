from pathlib import Path

import pytest
import yaml

from truemembers.config import (
    ConfigError,
    DataConfig,
    DomainsConfig,
    RunConfig,
    RunSection,
    SequenceSearchConfig,
    StructuralFilterConfig,
    StructureSearchConfig,
    load_config,
)


def make_sequence_db(tmp_path: Path) -> str:
    path = tmp_path / "uniprot.fasta"
    path.write_text(">sp|P00533|EGFR_HUMAN\nMKTL\n")
    return str(path)


def make_foldseek_db(tmp_path: Path) -> str:
    path = tmp_path / "afdb50"
    if not path.exists():
        path.write_text("")
    return str(path)


def make_hmms_folder(tmp_path: Path) -> str:
    path = tmp_path / "hmms"
    path.mkdir(exist_ok=True)
    return str(path)


def base_config(tmp_path: Path) -> dict:
    return {
        "run": {
            "target_uid": "Q5SK37",
            "target_region_start": 507,
            "target_region_end": 569,
            "output_dir": "runs",
            "convergence_mode": True,
            "min_growth": 0.2,
            "max_iterations": 10,
            "bins": 5,
            "hits_per_iteration": 10,
            "max_parallel_jobs": 4,
            "max_job_retries": 1,
            "cpu": 4,
        },
        "sequence_search": {
            "bit_scores": [0.2, 0.4, 0.6, 0.8],
            "iterations": 5,
        },
        "structure_search": {
            "mode": "local",
            "structure_db_path": make_foldseek_db(tmp_path),
            "tm_score_cutoff": 0.5,
            "max_hits": 1000,
        },
        "structural_filter": {
            "tool": "usalign",
            "tm_score_cutoff": 0.6,
        },
        "domains": {
            "hmms_folder_path": make_hmms_folder(tmp_path),
            "require_domains": False,
            "required_hmms": [],
        },
        "data": {
            "sequence_db_path": make_sequence_db(tmp_path),
            "sequence_db_index_dir": "",
            "structures_dir": "",
        },
    }


def write_config(tmp_path: Path, config: dict) -> Path:
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    return path


class TestLoadConfig:
    def test_loads_full_config(self, tmp_path):
        config = base_config(tmp_path)
        loaded = load_config(write_config(tmp_path, config))

        assert isinstance(loaded, RunConfig)
        assert loaded.run == RunSection(
            target_uid="Q5SK37",
            target_region_start=507,
            target_region_end=569,
            output_dir="runs",
            convergence_mode=True,
            min_growth=0.2,
            max_iterations=10,
            bins=5,
            hits_per_iteration=10,
            max_parallel_jobs=4,
            max_job_retries=1,
            cpu=4,
        )
        assert loaded.sequence_search == SequenceSearchConfig(
            bit_scores=[0.2, 0.4, 0.6, 0.8], iterations=5
        )
        assert loaded.structure_search == StructureSearchConfig(
            mode="local",
            structure_db_path=config["structure_search"]["structure_db_path"],
            tm_score_cutoff=0.5,
            max_hits=1000,
        )
        assert loaded.structural_filter == StructuralFilterConfig(
            tool="usalign", tm_score_cutoff=0.6
        )
        assert loaded.domains == DomainsConfig(
            hmms_folder_path=config["domains"]["hmms_folder_path"],
            require_domains=False,
            required_hmms=[],
        )
        assert loaded.data == DataConfig(
            sequence_db_path=config["data"]["sequence_db_path"],
            sequence_db_index_dir="",
            structures_dir="",
        )

    def test_optional_keys_fall_back_to_defaults(self, tmp_path):
        config = {
            "run": {
                "target_uid": "Q5SK37",
                "target_region_start": 507,
                "target_region_end": 569,
            },
            "sequence_search": {"bit_scores": [0.5]},
            "structure_search": {"mode": "webserver"},
            "domains": {"hmms_folder_path": ""},
            "data": {"sequence_db_path": make_sequence_db(tmp_path)},
        }
        loaded = load_config(write_config(tmp_path, config))

        assert loaded.run.output_dir == "runs"
        assert loaded.run.convergence_mode is True
        assert loaded.run.min_growth == 0.2
        assert loaded.run.max_iterations == 10
        assert loaded.run.bins == 5
        assert loaded.run.hits_per_iteration == 10
        assert loaded.run.max_parallel_jobs == 4
        assert loaded.run.max_job_retries == 1
        assert loaded.run.cpu == 4
        assert loaded.sequence_search.iterations == 5
        assert loaded.structure_search.tm_score_cutoff == 0.5
        assert loaded.structure_search.max_hits == 1000
        assert loaded.structural_filter.tool == "usalign"
        assert loaded.structural_filter.tm_score_cutoff == 0.6
        assert loaded.domains.require_domains is False
        assert loaded.domains.required_hmms == []
        assert loaded.data.sequence_db_index_dir == ""
        assert loaded.data.structures_dir == ""

    def test_bit_scores_preserve_int_and_float(self, tmp_path):
        config = base_config(tmp_path)
        config["sequence_search"]["bit_scores"] = [0.5, 10]
        loaded = load_config(write_config(tmp_path, config))

        assert loaded.sequence_search.bit_scores == [0.5, 10]
        assert isinstance(loaded.sequence_search.bit_scores[0], float)
        assert isinstance(loaded.sequence_search.bit_scores[1], int)

    def test_required_hmms_are_parsed(self, tmp_path):
        config = base_config(tmp_path)
        config["domains"]["require_domains"] = True
        config["domains"]["required_hmms"] = ["PF00069", "PF00512"]
        loaded = load_config(write_config(tmp_path, config))

        assert loaded.domains.require_domains is True
        assert loaded.domains.required_hmms == ["PF00069", "PF00512"]

    def test_to_dict_returns_config_snapshot(self, tmp_path):
        config = base_config(tmp_path)
        loaded = load_config(write_config(tmp_path, config))
        snapshot = loaded.to_dict()

        assert snapshot["run"]["target_uid"] == "Q5SK37"
        assert snapshot["run"]["target_region_start"] == 507
        assert snapshot["sequence_search"]["bit_scores"] == [0.2, 0.4, 0.6, 0.8]
        assert snapshot["structure_search"]["mode"] == "local"
        assert snapshot["domains"]["required_hmms"] == []
        assert snapshot["data"]["sequence_db_path"] == config["data"]["sequence_db_path"]

    def test_config_file_accepts_path_object(self, tmp_path):
        path = write_config(tmp_path, base_config(tmp_path))
        loaded = load_config(path)
        assert loaded.run.target_uid == "Q5SK37"


class TestMissingAndUnknownKeys:
    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(ConfigError, match="Config file not found"):
            load_config(tmp_path / "does_not_exist.yaml")

    def test_empty_file_raises(self, tmp_path):
        path = tmp_path / "empty.yaml"
        path.write_text("")
        with pytest.raises(ConfigError, match="empty"):
            load_config(path)

    def test_invalid_yaml_raises(self, tmp_path):
        path = tmp_path / "broken.yaml"
        path.write_text("run: [\n")
        with pytest.raises(ConfigError, match="Invalid YAML"):
            load_config(path)

    def test_root_must_be_mapping(self, tmp_path):
        path = tmp_path / "list.yaml"
        path.write_text("- a\n- b\n")
        with pytest.raises(ConfigError, match="mapping"):
            load_config(path)

    def test_unknown_top_level_key_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["extra"] = {}
        with pytest.raises(ConfigError, match="Unknown key"):
            load_config(write_config(tmp_path, config))

    def test_unknown_section_key_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["run"]["seed_uid"] = "P00533"
        with pytest.raises(ConfigError, match="Unknown key"):
            load_config(write_config(tmp_path, config))

    def test_section_must_be_mapping(self, tmp_path):
        config = base_config(tmp_path)
        config["run"] = "not-a-mapping"
        with pytest.raises(ConfigError, match="must be a mapping"):
            load_config(write_config(tmp_path, config))

    @pytest.mark.parametrize(
        "section,key",
        [
            ("run", "target_uid"),
            ("run", "target_region_start"),
            ("run", "target_region_end"),
            ("sequence_search", "bit_scores"),
            ("data", "sequence_db_path"),
        ],
    )
    def test_missing_required_key_raises(self, tmp_path, section, key):
        config = base_config(tmp_path)
        del config[section][key]
        with pytest.raises(ConfigError, match="Missing required key"):
            load_config(write_config(tmp_path, config))


class TestValueValidation:
    @pytest.mark.parametrize(
        "section,key,value,message",
        [
            ("run", "target_region_start", 0, "must be >= 1"),
            ("run", "target_region_end", 0, "must be >= 1"),
            ("run", "min_growth", -0.1, "must be >= 0"),
            ("run", "max_iterations", 0, "must be >= 1"),
            ("run", "bins", 0, "must be >= 1"),
            ("run", "hits_per_iteration", 0, "must be >= 1"),
            ("run", "max_parallel_jobs", 0, "must be >= 1"),
            ("run", "max_job_retries", -1, "must be >= 0"),
            ("run", "cpu", 0, "must be >= 1"),
            ("sequence_search", "iterations", 0, "must be >= 1"),
            ("structure_search", "max_hits", 0, "must be >= 1"),
            ("structure_search", "tm_score_cutoff", 1.5, "must be <= 1"),
            ("structure_search", "tm_score_cutoff", -0.5, "must be >= 0"),
            ("structural_filter", "tm_score_cutoff", 2, "must be <= 1"),
        ],
    )
    def test_out_of_range_values_raise(self, tmp_path, section, key, value, message):
        config = base_config(tmp_path)
        config[section][key] = value
        with pytest.raises(ConfigError, match=message):
            load_config(write_config(tmp_path, config))

    def test_region_end_before_start_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["run"]["target_region_end"] = 500
        with pytest.raises(ConfigError, match="target_region_end"):
            load_config(write_config(tmp_path, config))

    @pytest.mark.parametrize(
        "section,key,value",
        [
            ("run", "cpu", "four"),
            ("run", "cpu", 4.5),
            ("run", "max_iterations", True),
            ("run", "target_region_start", "507"),
        ],
    )
    def test_non_integer_values_raise(self, tmp_path, section, key, value):
        config = base_config(tmp_path)
        config[section][key] = value
        with pytest.raises(ConfigError, match="must be an integer"):
            load_config(write_config(tmp_path, config))

    @pytest.mark.parametrize(
        "section,key,value",
        [
            ("run", "convergence_mode", "yes"),
            ("domains", "require_domains", 1),
        ],
    )
    def test_non_boolean_values_raise(self, tmp_path, section, key, value):
        config = base_config(tmp_path)
        config[section][key] = value
        with pytest.raises(ConfigError, match="must be a boolean"):
            load_config(write_config(tmp_path, config))

    @pytest.mark.parametrize(
        "value",
        [[], "0.2", [0], [-0.1], ["high"], [True]],
    )
    def test_invalid_bit_scores_raise(self, tmp_path, value):
        config = base_config(tmp_path)
        config["sequence_search"]["bit_scores"] = value
        with pytest.raises(ConfigError, match="bit_scores"):
            load_config(write_config(tmp_path, config))

    @pytest.mark.parametrize(
        "value",
        ["PF00069", [1], [""], [None]],
    )
    def test_invalid_required_hmms_raise(self, tmp_path, value):
        config = base_config(tmp_path)
        config["domains"]["required_hmms"] = value
        with pytest.raises(ConfigError, match="required_hmms"):
            load_config(write_config(tmp_path, config))

    def test_non_string_target_uid_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["run"]["target_uid"] = 123
        with pytest.raises(ConfigError, match="must be a string"):
            load_config(write_config(tmp_path, config))

    def test_empty_target_uid_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["run"]["target_uid"] = ""
        with pytest.raises(ConfigError, match="non-empty"):
            load_config(write_config(tmp_path, config))


class TestModeAndToolValidation:
    def test_local_mode_requires_structure_db_path(self, tmp_path):
        config = base_config(tmp_path)
        config["structure_search"]["structure_db_path"] = ""
        with pytest.raises(ConfigError, match="structure_db_path"):
            load_config(write_config(tmp_path, config))

    def test_webserver_mode_allows_missing_structure_db_path(self, tmp_path):
        config = base_config(tmp_path)
        config["structure_search"]["mode"] = "webserver"
        config["structure_search"]["structure_db_path"] = ""
        loaded = load_config(write_config(tmp_path, config))
        assert loaded.structure_search.mode == "webserver"
        assert loaded.structure_search.structure_db_path == ""

    def test_unknown_structure_search_mode_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["structure_search"]["mode"] = "remote"
        with pytest.raises(ConfigError, match="mode"):
            load_config(write_config(tmp_path, config))

    def test_unknown_structural_filter_tool_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["structural_filter"]["tool"] = "both"
        with pytest.raises(ConfigError, match="tool"):
            load_config(write_config(tmp_path, config))


class TestPathValidation:
    def test_missing_sequence_db_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["data"]["sequence_db_path"] = str(tmp_path / "missing.fasta")
        with pytest.raises(ConfigError, match="sequence_db_path"):
            load_config(write_config(tmp_path, config))

    def test_sequence_db_must_be_file(self, tmp_path):
        config = base_config(tmp_path)
        config["data"]["sequence_db_path"] = str(tmp_path)
        with pytest.raises(ConfigError, match="sequence_db_path"):
            load_config(write_config(tmp_path, config))

    def test_missing_structure_db_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["structure_search"]["structure_db_path"] = str(tmp_path / "missing_db")
        with pytest.raises(ConfigError, match="structure_db_path"):
            load_config(write_config(tmp_path, config))

    def test_missing_hmms_folder_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["domains"]["hmms_folder_path"] = str(tmp_path / "missing_hmms")
        with pytest.raises(ConfigError, match="hmms_folder_path"):
            load_config(write_config(tmp_path, config))

    def test_require_domains_without_hmms_folder_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["domains"]["hmms_folder_path"] = ""
        config["domains"]["require_domains"] = True
        with pytest.raises(ConfigError, match="hmms_folder_path"):
            load_config(write_config(tmp_path, config))

    def test_empty_hmms_folder_disables_domains(self, tmp_path):
        config = base_config(tmp_path)
        config["domains"]["hmms_folder_path"] = ""
        loaded = load_config(write_config(tmp_path, config))
        assert loaded.domains.hmms_folder_path == ""

    def test_missing_structures_dir_raises(self, tmp_path):
        config = base_config(tmp_path)
        config["data"]["structures_dir"] = str(tmp_path / "missing_structures")
        with pytest.raises(ConfigError, match="structures_dir"):
            load_config(write_config(tmp_path, config))

    def test_existing_structures_dir_is_accepted(self, tmp_path):
        structures = tmp_path / "structures"
        structures.mkdir()
        config = base_config(tmp_path)
        config["data"]["structures_dir"] = str(structures)
        loaded = load_config(write_config(tmp_path, config))
        assert loaded.data.structures_dir == str(structures)

    def test_sequence_db_index_dir_may_not_exist_yet(self, tmp_path):
        config = base_config(tmp_path)
        config["data"]["sequence_db_index_dir"] = str(tmp_path / "future_index")
        loaded = load_config(write_config(tmp_path, config))
        assert loaded.data.sequence_db_index_dir == str(tmp_path / "future_index")

    def test_sequence_db_index_dir_must_be_directory_when_present(self, tmp_path):
        existing_file = tmp_path / "not_a_dir"
        existing_file.write_text("x")
        config = base_config(tmp_path)
        config["data"]["sequence_db_index_dir"] = str(existing_file)
        with pytest.raises(ConfigError, match="sequence_db_index_dir"):
            load_config(write_config(tmp_path, config))
