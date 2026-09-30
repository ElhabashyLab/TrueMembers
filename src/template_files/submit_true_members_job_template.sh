#!/bin/bash -l
#SBATCH --job-name=job_<seed_uid>_<iteration>
#SBATCH -o ./<run_dir>/job_<seed_uid>_<iteration>_%j.out
#SBATCH -e ./<run_dir>/job_<seed_uid>_<iteration>_%j.err
#SBATCH -D ./

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20000M
#SBATCH --time=7:00:00

cd <project_root>

module purge
module load evcouplings/0.2.1
module load python-waterboa/2024.06
eval "$(conda shell.bash hook)"
conda activate truemembers

python -u -m True_Members.run_job <job_config_path> <ev_couplings_config_template_path>