import logging
import json
import os
import subprocess
import time
from pathlib import Path

import yaml

from src.scripts.JobResultCollector import JobResultCollector
from src.scripts.JobStatus import JobStatus
from src.scripts.utils import extract_region_of_protein, run_download_tasks
from src.scripts.pipeline_db import (
    insert_member,
    insert_used_seed,
    mark_seed_failed,
    pick_unused_seeds,
)
from src.Hit import Hit, DetectionMethod
from src.Seed import Seed

logger = logging.getLogger(__name__)


def _parse_domains_json(raw: str) -> list[dict[str, tuple[int, int]]]:
    if not raw or raw == "[]":
        return []
    parsed = json.loads(raw)
    return [{k: tuple(v) for k, v in d.items()} for d in parsed]


class JobScheduler:
    
    def __init__(self, config, paths):
        self.start_config = config
        self.CONVERGENCE_RUN = self.start_config["global"]["convergence_mode"]
        self.MIN_GROWTH = self.start_config["global"]["min_growth"]
        self.MAX_ITERATIONS = self.start_config["global"]["max_iterations"]
        self.B = self.start_config["global"]["bins"]
        self.H = self.start_config["global"]["hits"]

        self.global_run_dir = paths["GLOBAL_RUN_DIR"]
        self.shared_structures_dir = paths["SHARED_STRUCTURES_DIR"]
        self.shared_sequences_dir = paths["SHARED_SEQUENCES_DIR"]
        self.results_dir = paths["RESULTS_DIR"]
        self.db_path = paths["DB_PATH"]

        self.job_configs_dir = paths["JOB_CONFIGS_DIR"]
        self.job_slurm_scripts_dir = paths["JOB_SLURM_SCRIPTS_DIR"]
        self.evc_config_template_path = paths["EVC_CONFIG_TEMPLATE"]
        self.submit_job_slurm_template_path = paths["SUBMIT_JOB_SLURM_TEMPLATE"]

        self.current_iteration = 1
        self.project_root = os.getcwd()
        self.jobs: dict[int, dict] = {} # {job_id: {"status": JobStatus, "seed": Seed}}
        self.job_results_collector = JobResultCollector(self.db_path)

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
        # temporary method: seeds are choosen arbitrarily from the database.
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
        job_config, job_config_path = self._build_job_config(hit)
        slurm_script_path = self._build_slurm_script(job_config, job_config_path)
        return Seed(hit=hit, job_config_path=job_config_path, slurm_script_path=slurm_script_path)


    def _build_job_config(self, hit: Hit) -> tuple[dict, Path]:
        root_dir = str(
            self.global_run_dir / "iterations" / "iteration_{self.current_iteration}" /
            f"seed_{hit.unique_id}_uid_{hit.uniprot_id}_{time.strftime('%Y-%m-%d_%H-%M-%S')}"
        )
        contents = {
            "seed": {
                "unique_id": hit.unique_id,
                "uniprot_id": hit.uniprot_id,
                "region": [hit.region_start, hit.region_end],
                "sequence_path": str(hit.sequence_path),
                "structure_path": str(hit.structure_path),
                "qtmscore": hit.qtmscore,
                "is_member": hit.is_member,
            },
            "pipeline": {
                "iteration": self.current_iteration,
                "root_dir": root_dir,
                "bit_scores": self.start_config["run"]["bit_scores"],
                "tm_score_cutoff": self.start_config["run"]["tm_score_cutoff"],
                "hmms_folder": self.start_config["run"]["hmms_folder_path"],
                "structural_filter_tool": self.start_config["run"].get("structural_filter_tool", "usalign"),
                "cpu": self.start_config["run"]["CPU"],
                "uniprot_fasta_path": self.start_config["run"].get("uniprot_fasta_path", ""),
                "uniprot_fasta_index_dir": self.start_config["run"].get("uniprot_fasta_index_dir", ""),
            },
            "shared": {
                "structures_dir": str(self.shared_structures_dir),
                "sequences_dir": str(self.shared_sequences_dir),
                "target_roi_structure": self.roi_of_initial_seed_structure,
            },
        }
        
        job_config_path = self.job_configs_dir / f"iteration_{self.current_iteration} / seed_{hit.unique_id}_uid_{hit.uniprot_id}.yaml"
        with open(job_config_path, "w") as f:
            yaml.safe_dump(contents, f, sort_keys=False)

        return contents, job_config_path


    def _build_slurm_script(self, job_config: dict, job_config_path: Path) -> Path:
        with open(self.submit_job_slurm_template_path) as f:
            slurm_template = f.read()
        
        slurm_script = slurm_template.replace("<iteration>", str(self.current_iteration))
        slurm_script = slurm_script.replace("<unique_id>", job_config["seed"]["unique_id"])
        slurm_script = slurm_script.replace("<uniprot_id>", job_config["seed"]["uniprot_id"])
        slurm_script = slurm_script.replace("<root_dir>", job_config["pipeline"]["root_dir"])
        slurm_script = slurm_script.replace("<project_root>", self.project_root)
        slurm_script = slurm_script.replace("<job_config_path>", str(job_config_path))
        slurm_script = slurm_script.replace("<ev_couplings_config_template_path>",
                                            str(self.evc_config_template_path))

        script_path = self.job_slurm_scripts_dir / f"iteration_{self.current_iteration} / seed_{job_config['seed']['unique_id']}_uid_{job_config['seed']['uniprot_id']}.sh"
        with open(script_path, "w") as f:
            f.write(slurm_script)

        return script_path
        

    def _start_job(self, seed: Seed, retried: bool = False) -> int:
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
            "seed": seed,
            "iteration": self.current_iteration,
            "retried": retried,
        }
        logger.info("Submitted job %d for seed %s (iteration %d).",
                     job_id, seed.unique_id, self.current_iteration)
        self._add_seed_to_used_seeds(seed.hit)
        return job_id


    def get_jobs(self, status: JobStatus) -> list[int]:
        return [job_id for job_id, info in self.jobs.items() if info["status"] == status]
    
    
    def _add_seed_to_used_seeds(self, hit: Hit):
        insert_used_seed(self.db_path, hit.unique_id, hit.uniprot_id, self.current_iteration)


    def _monitor_jobs(self):
        while True:
            self.update_job_status()
            self.job_results_collector.get_results_from_finished_jobs(self.jobs)
            self.restart_failed_jobs()
            
            pending = self.get_jobs(JobStatus.PENDING)
            running = self.get_jobs(JobStatus.RUNNING)
            active_jobs = pending + running
            
            if not active_jobs:
                logger.info("Iteration %d: all jobs finished.", self.current_iteration)
                break

            logger.info("Waiting for %d job(s) to finish...", len(active_jobs))
            time.sleep(30)

    
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
                    # Job not found -> treat as failed
                    self.jobs[job_id]["status"] = JobStatus.FAILED
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


    def restart_failed_jobs(self):
        failed_job_ids = self.get_jobs(JobStatus.FAILED)
        if not failed_job_ids:
            return

        for job_id in failed_job_ids:
            try:
                job_info = self.jobs[job_id]
                seed = job_info["seed"]

                if job_info.get("retried"):
                    logger.warning(
                        "Seed '%s' (uid %s, job %d) failed for the second time. Dropping job.",
                        seed.hit.uniprot_id, seed.hit.unique_id, job_id,
                    )
                    mark_seed_failed(self.db_path, seed.hit.unique_id)
                    job_info["status"] = JobStatus.HARVESTED
                    continue

                logger.info(
                    "Restarting job %d for seed '%s' (uid %s).",
                    job_id, seed.hit.uniprot_id, seed.hit.unique_id,
                )
                job_info["status"] = JobStatus.HARVESTED
                self._start_job(seed, retried=True)
            except Exception as e:
                logger.error("Failed to handle failed job %d: %s", job_id, e)


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
        logger.info("Pipeline ended. Results saved to %s.", self.db_path)


    def compute_convergence(self, absolute_floor: int, patience_window: int) -> tuple[bool, int]:
        import sqlite3
        with sqlite3.connect(str(self.db_path)) as conn:
            rows = conn.execute(
                "SELECT new_members FROM used_seeds "
                "WHERE failed = 0 AND new_members IS NOT NULL "
                "ORDER BY iteration"
            ).fetchall()
        finished = [r[0] for r in rows]
        if len(finished) < patience_window:
            below = sum(1 for v in finished if v < absolute_floor)
            return False, below
        window = finished[-patience_window:]
        below = sum(1 for v in window if v < absolute_floor)
        return all(v < absolute_floor for v in window), below
    