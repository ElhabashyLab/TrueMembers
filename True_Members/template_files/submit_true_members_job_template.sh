#!/bin/bash -l
#SBATCH --job-name=TRUE_MEMBERS_job_<seed_uid>_<iteration>         # Job name
#SBATCH -o ./<run_dir>/TRUE_MEMBERS_run_<seed_uid>_<iteration>_%j.out            # Output log file 
#SBATCH -e ./<run_dir>/TRUE_MEMBERS_run_<seed_uid>_<iteration>_%j.err            # Error log file
#SBATCH -D ./                                  # Initial working directory

#SBATCH --ntasks=1                             # Number of tasks (processes)
#SBATCH --cpus-per-task=10                   # Number of CPU cores for multithreading  
#SBATCH --mem=20000M                           # requested memory
#SBATCH --time=6:00:00                         # time limit

cd /ptmp/arowi/STAC_CbrA_Aggregation/true_members_pipeline_runs

module purge
module load python-waterboa/2024.06
module load evcouplings/0.2.1
eval "$(conda shell.bash hook)"
conda activate ba_thesis_new

python run_job.py <job_config_path> <ev_couplings_config_template_path>
