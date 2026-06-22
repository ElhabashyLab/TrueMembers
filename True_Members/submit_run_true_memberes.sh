#!/bin/bash -l
#SBATCH --job-name=true_members_pipeline_seed_Q5SK37
#SBATCH -o ./true_members_pipeline_seed_Q5SK37_%j.out
#SBATCH -e ./true_members_pipeline_seed_Q5SK37_%j.err
#SBATCH -D ./

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40000M
#SBATCH --time=24:00:00

# Change to the directory where sbatch was submitted (project root)
cd $SLURM_SUBMIT_DIR

module purge
module load evcouplings/0.2.1
module load python-waterboa/2024.06
eval "$(conda shell.bash hook)"
conda activate truemembers

truemembers run_true_members_config.yaml