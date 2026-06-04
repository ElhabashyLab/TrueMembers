import os
import subprocess
import yaml


def run_evcouplings(seed: list, bit_scores: list[float], cpus: int, config_template: str, root_dir: str) -> None:    
    config_file = set_up_config(seed, cpus, bit_scores, config_template, f"{root_dir}/data/evcouplings")        
    subprocess.run(["evcouplings",
                    config_file,
                    "--yolo"])
    
    try:
        result = subprocess.run(
            [
            "evcouplings", config_file, "--yolo"
            ],
            check=True,
            capture_output=True,
            text=True,)

    except subprocess.CalledProcessError as e:
        print("Foldseek failed!")
        print("Exit code:", e.returncode)
        print("STDOUT:")
        print(e.stdout)
        print("STDERR:")
        print(e.stderr)
        raise


def set_up_config(seed: list, cpus: int, bit_scores: list[float], config_template: str, evcouplings_dir: str) -> str:
    with open(config_template, "r") as f:
        config = yaml.safe_load(f)
        
    # set up prefix
    protein_uid = seed[0]
    r_start = seed[1]
    r_end = seed[2]
    uid_start_end = f"{protein_uid}_{r_start}-{r_end}"
    prefix = f"output/{uid_start_end}"
        
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
