import logging
import os
import subprocess

import yaml

from scripts.utils import read_yaml_config

logger = logging.getLogger(__name__)


def run_evcouplings(seed: list, bit_scores: list[float], cpus: int, config_template: str, root_dir: str) -> None:
    config_file = set_up_config(seed, cpus, bit_scores, config_template, f"{root_dir}/evcouplings")
    try:
        subprocess.run(
            [
                "evcouplings", config_file, "--yolo"
            ],
            check=True,
            capture_output=True,
            text=True,
        )

    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"evcouplings failed with error:\n{e.stderr}") from e


def set_up_config(seed: list, cpus: int, bit_scores: list[float], config_template: str, evcouplings_dir: str) -> str:
    config = read_yaml_config(config_template)

    protein_uid = seed[0]
    r_start = seed[1]
    r_end = seed[2]
    uid_start_end = f"{protein_uid}_{r_start}-{r_end}"
    prefix = f"{evcouplings_dir}/output/{uid_start_end}"

    config["global"]["cpu"] = cpus
    config["global"]["sequence_id"] = protein_uid
    config["global"]["region"] = [int(r_start), int(r_end)]
    config["global"]["prefix"] = prefix

    config["environment"]["cores"] = cpus
    config["environment"]["engine"] = "local"

    batch_config = {}
    for b in bit_scores:
        batch_key = f"_b{b}"
        batch_config[batch_key] = {"align": {"domain_threshold": b, "sequence_threshold": b}}

    if batch_config:
        config["batch"] = batch_config

    output_config_path = f"{evcouplings_dir}/config_{protein_uid}.yaml"
    with open(output_config_path, "w") as f:
        yaml.safe_dump(config, f)

    return output_config_path