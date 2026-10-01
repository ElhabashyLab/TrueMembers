import logging
import json
import os
import subprocess
import time

import pandas as pd
import yaml

from src.scripts.JobResultCollector import JobResultCollector
from src.scripts.JobStatus import JobStatus
from src.scripts.utils import copy_file, extract_region_of_protein, read_yaml_config, run_download_tasks
from src.scripts.pipeline_db import insert_member, insert_used_seed, pick_unused_seeds
from src.Hit import Hit, DetectionMethod
from src.Seed import Seed
from src.Target import Target

logger = logging.getLogger(__name__)


def _parse_domains_json(raw: str) -> list[dict[str, tuple[int, int]]]:
    if not raw or raw == "[]":
        return []
    parsed = json.loads(raw)
    return [{k: tuple(v) for k, v in d.items()} for d in parsed]


class JobScheduler:
    
    def __init__(self, config, paths):
        self.start_config = config
        self.CONVERGENCE_RUN = self.start_config["global"]["convergence_mode"] # if True, the pipeline will run until convergence
        self.MIN_GROWTH = self.start_config["global"]["min_growth"] # convergence criteria: minimum growth of Members in each iteration
        self.MAX_ITERATIONS = self.start_config["global"]["max_iterations"] # only used if not convergence_mode
        self.B = self.start_config["global"]["bins"] # number of bins which are used to seperate the sequence identity distribution
        self.H = self.start_config["global"]["hits"] # number of hits used as starting point in each iteration (chosen from bins)

        self.global_run_dir = paths["GLOBAL_RUN_DIR"]
        self.global_data_dir = paths["GLOBAL_DATA_DIR"]
        self.shared_structures_dir = paths["SHARED_STRUCTURES_DIR"]
        self.shared_sequences_dir = paths["SHARED_SEQUENCES_DIR"]
        self.results_dir = paths["RESULTS_DIR"]
    
        self.members_path = paths["MEMBERS_PATH"]
        self.used_seeds_path = paths["USED_SEEDS_PATH"]
        self.all_members_path = paths["ALL_MEMBERS_PATH"]
        self.sequence_identity_dir = paths["SEQUENCE_IDENTITY_DIR"]
        self.seed_to_members_sequence_identity_path = paths["SEED_TO_MEMBERS_SEQUENCE_IDENTITY_PATH"]
        self.db_path = paths["DB_PATH"]

        self.job_configs_dir = paths["JOB_CONFIGS_DIR"]
        self.job_slurm_scripts_dir = paths["JOB_SLURM_SCRIPTS_DIR"]
        self.job_config_template_path = paths["JOB_CONFIG_TEMPLATE"]
        self.evc_config_template_path = paths["EVC_CONFIG_TEMPLATE"]
        self.submit_job_slurm_template_path = paths["SUBMIT_JOB_SLURM_TEMPLATE"]

        self.current_iteration = 1
        self.project_root = os.getcwd()
        self._seeds_retried: set[str] = set()
        self.jobs: dict[int, dict] = {}
        self.job_results_collector = JobResultCollector(self.members_path,
                                                        self.used_seeds_path,
                                                        self.all_members_path,
                                                        self.sequence_identity_dir,
                                                        self.seed_to_members_sequence_identity_path,
                                                        self.shared_sequences_dir)

        self._setup_logging()


    def _setup_logging(self):
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.INFO)
        root_logger.handlers.clear()
        formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")

        pipeline_log = os.path.join(self.global_run_dir, "pipeline.log")
        fh = logging.FileHandler(pipeline_log)
        fh.setFormatter(formatter)
        root_logger.addHandler(fh)

        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        root_logger.addHandler(ch)
        
        
    def start(self):
        # The first iteration starts a single job with the target as seed
        target_seed = self._get_target_seed()
        initial_job_id = self._start_job(target_seed)
        self._monitor_jobs()
        
        # The main loop continues to start iterations with new seeds
        while True:
            seeds = self._pick_seeds()
            if not seeds:
                logger.info("No unused seeds remaining.")
                break
            
            ids_of_submitted_jobs = []
            for seed in seeds:
                job_id = self._start_job(seed)
                ids_of_submitted_jobs.append(job_id)
            logger.info("Iteration %d: submitted %d jobs, waiting...", self.current_iteration, len(ids_of_submitted_jobs))
            self._monitor_jobs()

            if self.CONVERGENCE_RUN and self._check_convergence():
                logger.info("Pipeline converged.")
                break

            if not self.CONVERGENCE_RUN and self.current_iteration > self.MAX_ITERATIONS:
                logger.info("Reached max iterations (%d).", self.MAX_ITERATIONS)
                break
            
            self.current_iteration += 1

        self.end_pipeline()
        
    
    def _get_target_seed(self) -> Seed:
        target_uid = self.start_config["global"]["target_uid"]
        region_start = self.start_config["global"]["target_region_start"]
        region_end = self.start_config["global"]["target_region_end"]

        run_download_tasks([target_uid], output_dir=str(self.shared_structures_dir), max_workers=1)
        structure_path = self.shared_structures_dir / f"{target_uid}.cif"
        self.roi_of_initial_seed_structure = extract_region_of_protein(
            str(structure_path), (region_start, region_end),
        )

        target = Hit(
            uniprot_id=target_uid,
            region_start=region_start,
            region_end=region_end,
            sequence_path=self.shared_sequences_dir / f"{target_uid}.fasta",
            structure_path=structure_path,
            detection_method=DetectionMethod.BOTH,
            qtmscore=1.0,
            is_member=True,
        )
        return self._create_seed(target)
    

    def _pick_seeds(self) -> list[Seed]:
        batch_size = self.B * self.H
        rows = pick_unused_seeds(self.db_path, batch_size)

        chosen_seeds: list[Seed] = []
        for row in rows:
            hit = self._row_to_hit(row)
            chosen_seeds.append(self._create_seed(hit))

        logger.info("Picked %d seeds for iteration %d.", len(chosen_seeds), self.current_iteration)
        return chosen_seeds

    def _row_to_hit(self, row: dict) -> Hit:
        return Hit(
            unique_id=row["unique_id"],
            uniprot_id=row["uniprot_id"],
            region_start=row["region_start"],
            region_end=row["region_end"],
            sequence_path=Path(row["sequence_path"]),
            structure_path=Path(row["structure_path"]),
            detection_method=DetectionMethod(row["detection_method"]),
            qtmscore=row["qtmscore"],
            is_member=bool(row["is_member"]),
            domains=_parse_domains_json(row.get("domains", "[]")),
        )


    def _create_seed(self, hit: Hit) -> Seed:
        job_config = read_yaml_config(str(self.job_config_template_path)).copy()
        job_config["run"]["unique_id"] = hit.unique_id
        job_config["run"]["seed_region_start"] = hit.region_start
        job_config["run"]["seed_region_end"] = hit.region_end
        job_config["global"]["iteration"] = self.current_iteration
        job_config["global"]["root_dir"] = str(
            self.global_run_dir / "runs" /
            f"job_{hit.unique_id}_iteration_{self.current_iteration}_{time.strftime('%Y-%m-%d_%H-%M-%S')}"
        )
        job_config["global"]["shared_structures_dir"] = str(self.shared_structures_dir)
        job_config["global"]["shared_sequences_dir"] = str(self.shared_sequences_dir)
        job_config["global"]["roi_structure_of_initial_seed"] = self.roi_of_initial_seed_structure
        job_config["run"]["bit_scores"] = self.start_config["run"]["bit_scores"]
        job_config["run"]["tm_score_cutoff"] = self.start_config["run"]["tm_score_cutoff"]
        job_config["run"]["hmms_folder_path"] = self.start_config["run"]["hmms_folder_path"]
        job_config["run"]["structural_filter_tool"] = self.start_config["run"].get("structural_filter_tool", "usalign")
        job_config["run"]["CPU"] = self.start_config["run"]["CPU"]
        job_config["run"]["uniprot_fasta_path"] = self.start_config["run"].get("uniprot_fasta_path", "")
        job_config["run"]["uniprot_fasta_index_dir"] = self.start_config["run"].get("uniprot_fasta_index_dir", "")

        seed_uid = hit.unique_id
        job_config_path = self.job_configs_dir / f"job_{seed_uid}_iteration_{self.current_iteration}.yaml"
        with open(job_config_path, "w") as f:
            yaml.safe_dump(job_config, f, sort_keys=False)

        slurm_script_path = self._create_slurm_script(str(job_config_path), job_config)

        self._add_seed_to_used_seeds(hit.unique_id, hit.uniprot_id)

        return Seed(hit=hit, job_config_path=job_config_path, slurm_script_path=slurm_script_path)


    def _start_job(self, seed: Seed) -> int:
        result = subprocess.run(
            ["sbatch", str(seed.slurm_script_path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        job_id = int(result.stdout.strip().split()[-1])
        self.jobs[job_id] = {
            "status": JobStatus.PENDING,
            "config_path": str(seed.job_config_path),
        }
        logger.info("Submitted job %d for seed %s (iteration %d).",
                     job_id, seed.unique_id, self.current_iteration)
        return job_id

    def _monitor_jobs(self):
        while True:
            self.update_job_status()
            self.job_results_collector.get_results_from_finished_jobs(self.jobs)
            self.remove_seeds_of_failed_jobs_from_used_seeds()

            pending_or_running = [
                jid for jid, j in self.jobs.items()
                if j["status"] in (JobStatus.PENDING, JobStatus.RUNNING)
            ]
            if not pending_or_running:
                logger.info("Iteration %d: all jobs finished.", self.current_iteration)
                break

            logger.info("Waiting for %d job(s) to finish...", len(pending_or_running))
            time.sleep(30)

    def _create_slurm_script(self, job_config_path: str, job_config: dict) -> Path:
        with open(self.submit_job_slurm_template_path) as f:
            slurm_template = f.read()

        seed_uid = job_config["run"]["seed_uid"]
        slurm_script = slurm_template.replace("<project_root>", self.project_root)
        slurm_script = slurm_script.replace("<run_dir>", job_config["global"]["root_dir"])
        slurm_script = slurm_script.replace("<seed_uid>", seed_uid)
        slurm_script = slurm_script.replace("<iteration>", str(self.current_iteration))
        slurm_script = slurm_script.replace("<job_config_path>", job_config_path)
        slurm_script = slurm_script.replace("<ev_couplings_config_template_path>",
                                            str(self.evc_config_template_path))

        script_path = self.job_slurm_scripts_dir / f"submit_job_{seed_uid}_iteration_{self.current_iteration}.sh"
        with open(script_path, "w") as f:
            f.write(slurm_script)

        return script_path

    def _add_seed_to_used_seeds(self, unique_id: str, uniprot_id: str):
        insert_used_seed(self.db_path, unique_id, uniprot_id, self.current_iteration)
    

    def download_roi_structure_of_initial_seed(self, start_config):
        target = start_config["global"]["target_uid"]
        region_start = start_config["global"]["target_region_start"]
        region_end = start_config["global"]["target_region_end"]
        run_download_tasks([target], output_dir=self.shared_structures_dir, max_workers=1)
        query_structure_full = f"{self.shared_structures_dir}/{target}.cif"
        query_structure_roi = extract_region_of_protein(query_structure_full, (region_start, region_end))
        return query_structure_roi


    def update_job_status(self) -> dict[int, dict]:
        for job_id in self.jobs.keys():
            try:
                result = subprocess.run(
                    [
                        "sacct",
                        "-j", str(job_id),
                        "-X",
                        "--format=State",
                        "--noheader",
                        "--parsable2"
                    ],
                    capture_output=True,
                    text=True,
                    check=True
                )
                state = result.stdout.strip().split()[0] if result.stdout.strip() else None
                if not state:
                    continue
                elif self.jobs[job_id]["status"] == JobStatus.HARVESTED:
                    continue
                elif state == "PENDING":
                    self.jobs[job_id]["status"] = JobStatus.PENDING
                elif state == "RUNNING":
                    self.jobs[job_id]["status"] = JobStatus.RUNNING
                elif state == "COMPLETED":
                    self.jobs[job_id]["status"] = JobStatus.COMPLETED
                elif state in ("FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL"):
                    self.jobs[job_id]["status"] = JobStatus.FAILED
                else:
                    logger.warning("Unknown job state '%s' for job %d. Keeping previous status %s.",
                                   state, job_id, self.jobs[job_id]['status'])

            except subprocess.CalledProcessError as e:
                logger.error("sacct error:\n%s", e.stderr)

        return self.jobs

    def remove_seeds_of_failed_jobs_from_used_seeds(self):
        failed_job_ids = [
            job_id for job_id in self.jobs.keys()
            if self.jobs[job_id]["status"] == JobStatus.FAILED
        ]
        if not failed_job_ids:
            return

        used_seeds_df = pd.read_csv(self.used_seeds_path)

        for job_id in failed_job_ids:
            try:
                job_config_path = self.jobs[job_id]["config_path"]
                job_config = read_yaml_config(job_config_path)
                failed_iteration = job_config["global"]["iteration"]
                seed_id_of_job = job_config["run"]["seed_uid"]

                if seed_id_of_job in self._seeds_retried:
                    logger.warning(
                        "Seed '%s' (iteration %d, job %d) failed for the second time. "
                        "Keeping it in used_seeds — no more retries.",
                        seed_id_of_job, failed_iteration, job_id
                    )
                    used_seeds_df.loc[used_seeds_df["iteration"] == failed_iteration, "new_true_members"] = 0
                    self.jobs[job_id]["status"] = JobStatus.HARVESTED
                    continue

                self._seeds_retried.add(seed_id_of_job)
                before = len(used_seeds_df)
                used_seeds_df = used_seeds_df[used_seeds_df["iteration"] != failed_iteration]
                removed = before - len(used_seeds_df)

                if removed:
                    logger.info("Freed seed '%s' (iteration %d) from used_seeds after job %d failed (retry 1/1).",
                                seed_id_of_job, failed_iteration, job_id)
                else:
                    logger.warning("No used_seeds row found for failed job %d (seed '%s', iteration %d).",
                                   job_id, seed_id_of_job, failed_iteration)
            except Exception as e:
                logger.error("Failed to handle failed job %d: %s", job_id, e)
            finally:
                if self.jobs[job_id]["status"] != JobStatus.HARVESTED:
                    self.jobs[job_id]["status"] = JobStatus.HARVESTED

        used_seeds_df.to_csv(self.used_seeds_path, index=False)

    def check_if_job_finished(self, job_id: int) -> bool:
        self.jobs = self.update_job_status()
        finished_jobs = (JobStatus.COMPLETED, JobStatus.HARVESTED, JobStatus.FAILED)
        return self.jobs[job_id]["status"] in finished_jobs

    def get_active_jobs_count(self) -> int:
        active_states = [JobStatus.RUNNING, JobStatus.PENDING]
        return len([job_id for job_id in self.jobs.keys() if self.jobs[job_id]["status"] in active_states])

    def end_pipeline(self):
        running_jobs_count = self.get_active_jobs_count()
        if running_jobs_count > 0:
            while True:
                self.jobs = self.update_job_status()
                running_jobs_count = self.get_active_jobs_count()

                if running_jobs_count == 0:
                    logger.info("All jobs finished. Ending pipeline.")
                    break
                else:
                    logger.info("Waiting for %d running jobs to finish before ending pipeline.", running_jobs_count)
                    time.sleep(180)

        self.jobs = self.update_job_status()
        self.job_results_collector.get_results_from_finished_jobs(self.jobs)
        logger.info("Saving final results to %s.", self.results_dir)
        copy_file(self.members_path, f"{self.results_dir}/True_Members_i_{self.current_iteration}.csv")
        copy_file(self.used_seeds_path, f"{self.results_dir}/Used_Seeds_i_{self.current_iteration}.csv")
        copy_file(self.all_members_path, f"{self.results_dir}/All_Members_i_{self.current_iteration}.csv")

    def compute_convergence(self, absolute_floor: int, patience_window: int) -> tuple[bool, int]:
        used_seeds_df = pd.read_csv(self.used_seeds_path)
        finished_seeds_df = used_seeds_df.dropna(subset=["new_true_members"]).sort_values(by="iteration").copy()

        if len(finished_seeds_df) < patience_window:
            jobs_in_window = finished_seeds_df
            jobs_in_window["below_floor"] = jobs_in_window["new_true_members"] < absolute_floor
            number_of_jobs_below_floor = int(jobs_in_window["below_floor"].sum())
            return False, number_of_jobs_below_floor
        else:
            jobs_in_window = finished_seeds_df.tail(patience_window).copy()
            jobs_in_window["below_floor"] = jobs_in_window["new_true_members"] < absolute_floor
            number_of_jobs_below_floor = int(jobs_in_window["below_floor"].sum())
            converged = bool(jobs_in_window["below_floor"].all())
            return converged, number_of_jobs_below_floor

    def start_job(self) -> int:
        job_config_path, job_config = self.create_job_config()
        slurm_script_path = self.create_slurm_script(job_config_path, job_config)
        self.add_current_seed_to_used_seeds(job_config)
        try:
            result = subprocess.run(
                ["sbatch", slurm_script_path],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True)

            job_id = int(result.stdout.strip().split()[-1])
            self.jobs[job_id] = {"status": JobStatus.PENDING, "config_path": job_config_path}
            logger.info("Submitted job %d for seed %s (iteration %d).", job_id,
                        job_config["run"]["seed_uid"], self.current_iteration)
            return job_id
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to submit job {slurm_script_path} with error:\n{e.stderr}") from e

    def add_current_seed_to_used_seeds(self, job_config):
        seed_row = pd.DataFrame([{
            "id": job_config["run"]["seed_uid"],
            "region_start": job_config["run"]["seed_region_start"],
            "region_end": job_config["run"]["seed_region_end"],
            "iteration": self.current_iteration,
            "new_true_members": pd.NA
        }])
        used_seeds_df = pd.read_csv(self.used_seeds_path)
        used_seeds_df = pd.concat([used_seeds_df, seed_row], ignore_index=True)
        used_seeds_df.to_csv(self.used_seeds_path, index=False)

    def create_slurm_script(self, job_config_path, job_config):
        with open(self.submit_job_slurm_template_path, "r") as f:
            slurm_template = f.read()

        seed_uid = job_config["run"]["seed_uid"]
        slurm_script = slurm_template.replace("<project_root>", self.project_root)
        slurm_script = slurm_script.replace("<run_dir>", job_config["global"]["root_dir"])
        slurm_script = slurm_script.replace("<seed_uid>", seed_uid)
        slurm_script = slurm_script.replace("<iteration>", str(self.current_iteration))
        slurm_script = slurm_script.replace("<job_config_path>", job_config_path)
        slurm_script = slurm_script.replace("<ev_couplings_config_template_path>", self.evc_config_template_path)

        job_slurm_script_path = f"{self.job_slurm_scripts_dir}/submit_job_{seed_uid}_iteration_{self.current_iteration}.sh"
        with open(job_slurm_script_path, "w") as f:
            f.write(slurm_script)

        return job_slurm_script_path


    def create_job_config(self):
        job_config = read_yaml_config(self.job_config_template_path).copy()

        if self.current_iteration == 1:
            new_seed = self.start_config["global"]["target_uid"]
            job_config["run"]["seed_uid"] = new_seed
            job_config["run"]["seed_region_start"] = self.start_config["global"]["target_region_start"]
            job_config["run"]["seed_region_end"] = self.start_config["global"]["target_region_end"]

        else:
            top_seed_row = self.get_highest_scoring_seed()
            new_seed = top_seed_row["id"].iloc[0]
            job_config["run"]["seed_uid"] = new_seed
            job_config["run"]["seed_region_start"] = int(top_seed_row["region_start"].iloc[0])
            job_config["run"]["seed_region_end"] = int(top_seed_row["region_end"].iloc[0])

        job_config["global"]["roi_structure_of_initial_seed"] = self.roi_of_initial_seed_structure
        job_config["global"]["iteration"] = self.current_iteration
        job_config["global"]["root_dir"] = f"{self.global_run_dir}/runs/job_{new_seed}_iteration_{self.current_iteration}_{time.strftime('%Y-%m-%d_%H-%M-%S')}"
        job_config["global"]["shared_structures_dir"] = self.shared_structures_dir
        job_config["global"]["shared_sequences_dir"] = self.shared_sequences_dir
        job_config["run"]["bit_scores"] = self.start_config["run"]["bit_scores"]
        job_config["run"]["tm_score_cutoff"] = self.start_config["run"]["tm_score_cutoff"]
        job_config["run"]["hmms_folder_path"] = self.start_config["run"]["hmms_folder_path"]
        job_config["run"]["structural_filter_tool"] = self.start_config["run"].get("structural_filter_tool", "usalign")
        job_config["run"]["CPU"] = self.start_config["run"]["CPU"]
        job_config["run"]["uniprot_fasta_path"] = self.start_config["run"].get("uniprot_fasta_path", "")
        job_config["run"]["uniprot_fasta_index_dir"] = self.start_config["run"].get("uniprot_fasta_index_dir", "")

        job_config_path = f"{self.job_configs_dir}/job_{new_seed}_iteration_{self.current_iteration}.yaml"
        with open(job_config_path, "w") as f:
            yaml.safe_dump(job_config, f, sort_keys=False)

        return job_config_path, job_config


    def get_highest_scoring_seed(self) -> pd.DataFrame:
        seeds_df = pd.read_csv(self.members_path)
        used_seeds_df = pd.read_csv(self.used_seeds_path)
        unused_seeds_df = seeds_df[~seeds_df["id"].isin(used_seeds_df["id"])]
        unused_seeds_df = unused_seeds_df.sort_values(by="qtmscore", ascending=False)
        top_seed_row = unused_seeds_df.iloc[[0]]
        return top_seed_row
    