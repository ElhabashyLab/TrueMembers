import logging
import os
import subprocess
import time

import pandas as pd
import yaml

from True_Members.scripts.JobResultCollector import JobResultCollector
from True_Members.scripts.JobStatus import JobStatus
from True_Members.scripts.utils import copy_file, extract_region_of_protein, read_yaml_config, run_download_tasks

logger = logging.getLogger(__name__)


class JobScheduler:
    def __init__(self,
                 config,
                 global_run_dir,
                 global_data_dir,
                 shared_structures_dir,
                 shared_sequences_dir,
                 seed_list_path,
                 used_seeds_path,
                 all_members_path,
                 results_dir,
                 job_configs_dir,
                 job_slurm_scripts_dir,
                 job_config_template_path,
                 evc_config_template_path,
                 submit_job_slurm_template_path):

        self.global_run_dir = global_run_dir
        self.shared_structures_dir = shared_structures_dir
        self.shared_sequences_dir = shared_sequences_dir
        self.start_config = config
        self.MAX_ITERATIONS = self.start_config["global"]["iterations"]
        self.MAX_PARALLEL_JOBS = self.start_config["global"]["max_parallel_jobs"]
        self.CONVERGENCE_RUN = self.start_config["global"]["convergence_run"]
        self.ABSOLUTE_FLOOR = self.start_config["global"]["absolute_floor"]
        self.PATIENCE_WINDOW = self.start_config["global"]["patience_window"]

        self.job_configs_dir = job_configs_dir
        self.job_slurm_scripts_dir = job_slurm_scripts_dir
        self.global_data_dir = global_data_dir
        self.seed_list_path = seed_list_path
        self.used_seeds_path = used_seeds_path
        self.all_members_path = all_members_path
        self.results_dir = results_dir

        self.job_config_template_path = job_config_template_path
        self.evc_config_template_path = evc_config_template_path
        self.submit_job_slurm_template_path = submit_job_slurm_template_path

        self.current_iteration = 1
        self.project_root = os.getcwd()
        self.jobs: dict[int, dict] = {}
        self.job_results_collector = JobResultCollector(self.seed_list_path, self.used_seeds_path, self.all_members_path)

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
        while True:
            self.jobs = self.update_job_status()
            self.job_results_collector.get_results_from_finished_jobs(self.jobs)
            self.remove_seeds_of_failed_jobs_from_used_seeds()
            active_jobs_count = self.get_active_jobs_count()
            if self.current_iteration != 1:
                seed_list = pd.read_csv(f"{self.seed_list_path}")
                used_seeds = pd.read_csv(f"{self.used_seeds_path}")
                available_seeds = seed_list[~seed_list["id"].isin(used_seeds["id"])]
                if not self.CONVERGENCE_RUN and self.current_iteration > self.MAX_ITERATIONS:
                    logger.info("Reached maximum number of iterations (%d). Stopping pipeline.", self.MAX_ITERATIONS)
                    self.end_pipeline()
                    break
                elif available_seeds.empty:
                    logger.info("No available seeds to run.")
                    if active_jobs_count > 0:
                        logger.info("Waiting for running jobs to finish. Currently running %d jobs.", active_jobs_count)
                        time.sleep(180)
                        continue
                    else:
                        logger.info("No running jobs. Stopping pipeline.")
                        self.end_pipeline()
                        break
                elif self.MAX_PARALLEL_JOBS <= active_jobs_count:
                    logger.info("Currently running %d/%d jobs. Waiting for jobs to finish.",
                                active_jobs_count, self.MAX_PARALLEL_JOBS)
                    time.sleep(180)
                    continue
                else:
                    if self.CONVERGENCE_RUN:
                        converged, number_below_floor = self.compute_convergence(self.ABSOLUTE_FLOOR, self.PATIENCE_WINDOW)
                        if converged:
                            logger.info("Number of True Members converged with %d consecutive jobs below absolute floor of %d. Stopping pipeline.",
                                        number_below_floor, self.ABSOLUTE_FLOOR)
                            self.end_pipeline()
                            break
                        else:
                            logger.info("Number of True Members did not converge with %d consecutive jobs below absolute floor of %d. Continuing pipeline.",
                                        number_below_floor, self.ABSOLUTE_FLOOR)
                            logger.info("Starting iteration %d until convergence.", self.current_iteration)
                            self.start_job()
                            self.current_iteration += 1
                    else:
                        logger.info("Starting iteration %d/%d.", self.current_iteration, self.MAX_ITERATIONS)
                        self.start_job()
                        self.current_iteration += 1
            else:
                logger.info("Starting initial iteration %d of %d", self.current_iteration, self.MAX_ITERATIONS)
                self.roi_of_initial_seed_structure = self.download_roi_structure_of_initial_seed(self.start_config)
                job_id = self.start_job()
                while True:
                    logger.info("Waiting for initial job to finish ...")
                    if self.check_if_job_finished(job_id):
                        logger.info("First iteration finished.")
                        self.current_iteration += 1
                        break
                    else:
                        time.sleep(180)

    def download_roi_structure_of_initial_seed(self, start_config):
        initial_seed = start_config["run"]["seed_uid"]
        region_start = start_config["run"]["seed_region_start"]
        region_end = start_config["run"]["seed_region_end"]
        run_download_tasks([initial_seed], output_dir=self.shared_structures_dir, max_workers=1)
        query_structure_full = f"{self.shared_structures_dir}/{initial_seed}.cif"
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

                before = len(used_seeds_df)
                used_seeds_df = used_seeds_df[used_seeds_df["iteration"] != failed_iteration]
                removed = before - len(used_seeds_df)

                if removed:
                    logger.info("Freed seed '%s' (iteration %d) from used_seeds after job %d failed.",
                                seed_id_of_job, failed_iteration, job_id)
                else:
                    logger.warning("No used_seeds row found for failed job %d (seed '%s', iteration %d).",
                                   job_id, seed_id_of_job, failed_iteration)
            except Exception as e:
                logger.error("Failed to free seed for failed job %d: %s", job_id, e)
            finally:
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
        copy_file(self.seed_list_path, f"{self.results_dir}/final_True_Members_i_{self.current_iteration}.csv")
        copy_file(self.used_seeds_path, f"{self.results_dir}/used_seeds_iteration_i_{self.current_iteration}.csv")
        copy_file(self.all_members_path, f"{self.results_dir}/all_members_iteration_i_{self.current_iteration}.csv")

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
            new_seed = self.start_config["run"]["seed_uid"]
            job_config["run"]["seed_uid"] = new_seed
            job_config["run"]["seed_region_start"] = self.start_config["run"]["seed_region_start"]
            job_config["run"]["seed_region_end"] = self.start_config["run"]["seed_region_end"]

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

        job_config_path = f"{self.job_configs_dir}/job_{new_seed}_iteration_{self.current_iteration}.yaml"
        with open(job_config_path, "w") as f:
            yaml.safe_dump(job_config, f, sort_keys=False)

        return job_config_path, job_config

    def get_highest_scoring_seed(self) -> pd.DataFrame:
        seeds_df = pd.read_csv(self.seed_list_path)
        used_seeds_df = pd.read_csv(self.used_seeds_path)
        available_seeds_df = seeds_df[~seeds_df["id"].isin(used_seeds_df["id"])]
        available_seeds_df = available_seeds_df.sort_values(by="qtmscore", ascending=False)
        top_seed_row = available_seeds_df.iloc[[0]]
        return top_seed_row
    