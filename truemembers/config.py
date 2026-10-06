from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field
from typing import Any

import yaml


class ConfigError(ValueError):
    """Raised when a configuration file is missing, malformed, or invalid."""


@dataclass(frozen=True)
class RunSection:
    target_uid: str
    target_region_start: int
    target_region_end: int
    output_dir: str = "runs"
    convergence_mode: bool = True
    min_growth: float = 0.2
    max_iterations: int = 10
    bins: int = 5
    hits_per_iteration: int = 10
    max_parallel_jobs: int = 4
    max_job_retries: int = 1
    cpu: int = 4


@dataclass(frozen=True)
class SequenceSearchConfig:
    bit_scores: list[int | float]
    iterations: int = 5


@dataclass(frozen=True)
class StructureSearchConfig:
    mode: str = "local"
    structure_db_path: str = ""
    tm_score_cutoff: float = 0.5
    max_hits: int = 1000


@dataclass(frozen=True)
class StructuralFilterConfig:
    tool: str = "usalign"
    tm_score_cutoff: float = 0.6


@dataclass(frozen=True)
class DomainsConfig:
    hmms_folder_path: str = "hmm_domain_filters"
    require_domains: bool = False
    required_hmms: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class DataConfig:
    sequence_db_path: str = ""
    sequence_db_index_dir: str = ""
    structures_dir: str = ""


@dataclass(frozen=True)
class RunConfig:
    run: RunSection
    sequence_search: SequenceSearchConfig
    structure_search: StructureSearchConfig
    structural_filter: StructuralFilterConfig
    domains: DomainsConfig
    data: DataConfig

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


_SECTIONS = {"run", "sequence_search", "structure_search", "structural_filter", "domains", "data"}

_RUN_KEYS = {
    "target_uid",
    "target_region_start",
    "target_region_end",
    "output_dir",
    "convergence_mode",
    "min_growth",
    "max_iterations",
    "bins",
    "hits_per_iteration",
    "max_parallel_jobs",
    "max_job_retries",
    "cpu",
}
_SEQUENCE_SEARCH_KEYS = {"bit_scores", "iterations"}
_STRUCTURE_SEARCH_KEYS = {"mode", "structure_db_path", "tm_score_cutoff", "max_hits"}
_STRUCTURAL_FILTER_KEYS = {"tool", "tm_score_cutoff"}
_DOMAINS_KEYS = {"hmms_folder_path", "require_domains", "required_hmms"}
_DATA_KEYS = {"sequence_db_path", "sequence_db_index_dir", "structures_dir"}

_STRUCTURE_SEARCH_MODES = ("local", "webserver")
_STRUCTURAL_FILTER_TOOLS = ("usalign", "foldseek")


def load_config(path: str | os.PathLike) -> RunConfig:
    raw = _read_config_file(path)
    _check_keys(raw, _SECTIONS, "config")
    return RunConfig(
        run=_load_run_section(raw),
        sequence_search=_load_sequence_search_section(raw),
        structure_search=_load_structure_search_section(raw),
        structural_filter=_load_structural_filter_section(raw),
        domains=_load_domains_section(raw),
        data=_load_data_section(raw),
    )


def _read_config_file(path: str | os.PathLike) -> dict:
    if not os.path.isfile(path):
        raise ConfigError(f"Config file not found: {path}")
    try:
        with open(path, "r") as file:
            raw = yaml.safe_load(file)
    except yaml.YAMLError as error:
        raise ConfigError(f"Invalid YAML in {path}: {error}") from error
    if raw is None:
        raise ConfigError(f"Config file is empty: {path}")
    if not isinstance(raw, dict):
        raise ConfigError("Config root must be a mapping of sections")
    return raw


def _load_run_section(raw: dict) -> RunSection:
    section = _section(raw, "run")
    _check_keys(section, _RUN_KEYS, "run")

    region_start = _get_int(section, "target_region_start", "run", required=True, minimum=1)
    region_end = _get_int(section, "target_region_end", "run", required=True, minimum=1)
    if region_end < region_start:
        raise ConfigError(
            "'run.target_region_end' must be >= 'run.target_region_start' "
            f"({region_end} < {region_start})"
        )

    return RunSection(
        target_uid=_get_str(section, "target_uid", "run", required=True),
        target_region_start=region_start,
        target_region_end=region_end,
        output_dir=_get_str(section, "output_dir", "run", default="runs"),
        convergence_mode=_get_bool(section, "convergence_mode", "run", default=True),
        min_growth=_get_float(section, "min_growth", "run", default=0.2, minimum=0.0),
        max_iterations=_get_int(section, "max_iterations", "run", default=10, minimum=1),
        bins=_get_int(section, "bins", "run", default=5, minimum=1),
        hits_per_iteration=_get_int(section, "hits_per_iteration", "run", default=10, minimum=1),
        max_parallel_jobs=_get_int(section, "max_parallel_jobs", "run", default=4, minimum=1),
        max_job_retries=_get_int(section, "max_job_retries", "run", default=1, minimum=0),
        cpu=_get_int(section, "cpu", "run", default=4, minimum=1),
    )


def _load_sequence_search_section(raw: dict) -> SequenceSearchConfig:
    section = _section(raw, "sequence_search")
    _check_keys(section, _SEQUENCE_SEARCH_KEYS, "sequence_search")
    return SequenceSearchConfig(
        bit_scores=_get_bit_scores(section, "sequence_search"),
        iterations=_get_int(section, "iterations", "sequence_search", default=5, minimum=1),
    )


def _load_structure_search_section(raw: dict) -> StructureSearchConfig:
    section = _section(raw, "structure_search")
    _check_keys(section, _STRUCTURE_SEARCH_KEYS, "structure_search")

    mode = _get_str(section, "mode", "structure_search", default="local")
    if mode not in _STRUCTURE_SEARCH_MODES:
        raise ConfigError(
            "'structure_search.mode' must be one of: "
            f"{', '.join(_STRUCTURE_SEARCH_MODES)} (got {mode!r})"
        )

    structure_db_path = _get_str(
        section, "structure_db_path", "structure_search", default="", allow_empty=True
    )
    if mode == "local" and not structure_db_path:
        raise ConfigError(
            "'structure_search.structure_db_path' is required when 'structure_search.mode' is 'local'"
        )
    if structure_db_path:
        _check_exists(structure_db_path, "structure_search.structure_db_path")

    return StructureSearchConfig(
        mode=mode,
        structure_db_path=structure_db_path,
        tm_score_cutoff=_get_float(
            section, "tm_score_cutoff", "structure_search", default=0.5, minimum=0.0, maximum=1.0
        ),
        max_hits=_get_int(section, "max_hits", "structure_search", default=1000, minimum=1),
    )


def _load_structural_filter_section(raw: dict) -> StructuralFilterConfig:
    section = _section(raw, "structural_filter")
    _check_keys(section, _STRUCTURAL_FILTER_KEYS, "structural_filter")

    tool = _get_str(section, "tool", "structural_filter", default="usalign")
    if tool not in _STRUCTURAL_FILTER_TOOLS:
        raise ConfigError(
            "'structural_filter.tool' must be one of: "
            f"{', '.join(_STRUCTURAL_FILTER_TOOLS)} (got {tool!r})"
        )

    return StructuralFilterConfig(
        tool=tool,
        tm_score_cutoff=_get_float(
            section, "tm_score_cutoff", "structural_filter", default=0.6, minimum=0.0, maximum=1.0
        ),
    )


def _load_domains_section(raw: dict) -> DomainsConfig:
    section = _section(raw, "domains")
    _check_keys(section, _DOMAINS_KEYS, "domains")

    require_domains = _get_bool(section, "require_domains", "domains", default=False)
    hmms_folder_path = _get_str(
        section, "hmms_folder_path", "domains", default="hmm_domain_filters", allow_empty=True
    )
    if not hmms_folder_path and require_domains:
        raise ConfigError(
            "'domains.hmms_folder_path' is required when 'domains.require_domains' is true"
        )
    if hmms_folder_path:
        _check_dir(hmms_folder_path, "domains.hmms_folder_path")

    return DomainsConfig(
        hmms_folder_path=hmms_folder_path,
        require_domains=require_domains,
        required_hmms=_get_str_list(section, "required_hmms", "domains", default=[]),
    )


def _load_data_section(raw: dict) -> DataConfig:
    section = _section(raw, "data")
    _check_keys(section, _DATA_KEYS, "data")

    sequence_db_path = _get_str(section, "sequence_db_path", "data", required=True)
    _check_file(sequence_db_path, "data.sequence_db_path")

    sequence_db_index_dir = _get_str(
        section, "sequence_db_index_dir", "data", default="", allow_empty=True
    )
    if sequence_db_index_dir:
        _check_dir_if_exists(sequence_db_index_dir, "data.sequence_db_index_dir")

    structures_dir = _get_str(section, "structures_dir", "data", default="", allow_empty=True)
    if structures_dir:
        _check_dir(structures_dir, "data.structures_dir")

    return DataConfig(
        sequence_db_path=sequence_db_path,
        sequence_db_index_dir=sequence_db_index_dir,
        structures_dir=structures_dir,
    )


def _section(raw: dict, name: str) -> dict:
    value = raw.get(name, {})
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ConfigError(f"Section '{name}' must be a mapping")
    return value


def _check_keys(section: dict, allowed: set[str], where: str) -> None:
    unknown = sorted(set(section) - allowed)
    if unknown:
        raise ConfigError(f"Unknown key(s) in '{where}': {', '.join(unknown)}")


def _get_str(
    section: dict,
    key: str,
    where: str,
    default: str | None = None,
    required: bool = False,
    allow_empty: bool = False,
) -> str:
    if key not in section:
        if required:
            raise ConfigError(f"Missing required key '{where}.{key}'")
        return default
    value = section[key]
    if not isinstance(value, str):
        raise ConfigError(f"'{where}.{key}' must be a string")
    if value == "":
        if allow_empty:
            return value
        raise ConfigError(f"'{where}.{key}' must be a non-empty string")
    if not value.strip():
        raise ConfigError(f"'{where}.{key}' must be a non-empty string")
    return value


def _get_bool(section: dict, key: str, where: str, default: bool) -> bool:
    value = section.get(key, default)
    if not isinstance(value, bool):
        raise ConfigError(f"'{where}.{key}' must be a boolean")
    return value


def _get_int(
    section: dict,
    key: str,
    where: str,
    default: int | None = None,
    required: bool = False,
    minimum: int | None = None,
) -> int:
    if key not in section:
        if required:
            raise ConfigError(f"Missing required key '{where}.{key}'")
        return default
    value = section[key]
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigError(f"'{where}.{key}' must be an integer")
    if minimum is not None and value < minimum:
        raise ConfigError(f"'{where}.{key}' must be >= {minimum} (got {value})")
    return value


def _get_float(
    section: dict,
    key: str,
    where: str,
    default: float | None = None,
    minimum: float | None = None,
    maximum: float | None = None,
) -> float:
    value = section.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"'{where}.{key}' must be a number")
    if minimum is not None and value < minimum:
        raise ConfigError(f"'{where}.{key}' must be >= {minimum} (got {value})")
    if maximum is not None and value > maximum:
        raise ConfigError(f"'{where}.{key}' must be <= {maximum} (got {value})")
    return float(value)


def _get_str_list(section: dict, key: str, where: str, default: list[str]) -> list[str]:
    value = section.get(key, default)
    if not isinstance(value, list):
        raise ConfigError(f"'{where}.{key}' must be a list of strings")
    for index, entry in enumerate(value):
        if not isinstance(entry, str) or not entry.strip():
            raise ConfigError(f"'{where}.{key}[{index}]' must be a non-empty string")
    return list(value)


def _get_bit_scores(section: dict, where: str) -> list[int | float]:
    if "bit_scores" not in section:
        raise ConfigError(f"Missing required key '{where}.bit_scores'")
    value = section["bit_scores"]
    if not isinstance(value, list) or not value:
        raise ConfigError(f"'{where}.bit_scores' must be a non-empty list of numbers")
    for index, entry in enumerate(value):
        if isinstance(entry, bool) or not isinstance(entry, (int, float)):
            raise ConfigError(f"'{where}.bit_scores[{index}]' must be a number")
        if entry <= 0:
            raise ConfigError(f"'{where}.bit_scores[{index}]' must be > 0 (got {entry})")
    return list(value)


def _check_exists(path: str, where: str) -> None:
    if not os.path.exists(path):
        raise ConfigError(f"'{where}' does not exist: {path}")


def _check_file(path: str, where: str) -> None:
    if not os.path.isfile(path):
        raise ConfigError(f"'{where}' must point to an existing file: {path}")


def _check_dir(path: str, where: str) -> None:
    if not os.path.isdir(path):
        raise ConfigError(f"'{where}' must point to an existing directory: {path}")


def _check_dir_if_exists(path: str, where: str) -> None:
    if os.path.exists(path) and not os.path.isdir(path):
        raise ConfigError(f"'{where}' must be a directory: {path}")
