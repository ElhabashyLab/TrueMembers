#!/bin/bash -l
#SBATCH --job-name=true_members_pipeline_seed_Q5SK37         # Job name
#SBATCH -o ./true_members_pipeline_seed_Q5SK37_%j.out            # Output log file 
#SBATCH -e ./true_members_pipeline_seed_Q5SK37_%j.err            # Error log file
#SBATCH -D ./                                  # Initial working directory

#SBATCH --ntasks=1                             # Number of tasks (processes)
#SBATCH --cpus-per-task=5                   # Number of CPU cores for multithreading  
#SBATCH --mem=40000M                           # requested memory
#SBATCH --time=24:00:00                         # time limit

cd /ptmp/arowi/STAC_CbrA_Aggregation/true_members_pipeline_runs

module purge
module load python-waterboa/2024.06
eval "$(conda shell.bash hook)"
conda activate ba_thesis_new

python -u run_true_members.py run_true_members_config.yaml
