#!/bin/sh

#SBATCH --time=2-00:00:00
#SBATCH --job-name=fbr_edit
#SBATCH --error=sens_fbr_edit.log
#SBATCH --output=sens_fbr_edit.slurm.log
#SBATCH --partition=short

python simulation_custom_Tprofiles_sensitivity.py

