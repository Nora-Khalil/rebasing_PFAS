#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=17:00:00
#SBATCH --job-name=all_fixing_reg
#SBATCH --error=fixing_reg_error_%a.slurm.log
#SBATCH --output=fixing_reg_%a.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --array=12,17
#SBATCH --partition=short

index=$((SLURM_ARRAY_TASK_ID - 1))


python retraining_all_PFAS_fams.py $index