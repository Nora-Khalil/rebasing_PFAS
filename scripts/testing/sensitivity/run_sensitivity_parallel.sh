#!/bin/bash
#SBATCH --job-name=sensitivity
#SBATCH --output=diagram_all.slurm.%x.log
#SBATCH --error=diagram_fc_all.slurm.%x.log
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1

#list_of_blends=(C2H5F_CH2F2 C2H5F_CH2FCH2F C2H5F_CH2FCHF2 C2H5F_CH3CHF2 C2H5F_CH3F CH2F2_CH3CF3 CH2F2_CH3F CH2FCH2F_CH2F2 CH3CF3_C2H5F CH3CF3_CH2FCH2F CH3CF3_CH2FCHF2 CH3CF3_CH3CHF2 CH3CHF2_CH2F2 CH3CHF2_CH2FCHF2 CH3F_CH3CF3)
list_of_blends=(CH3F)

index=$SLURM_ARRAY_TASK_ID-1

blend_name="${list_of_blends[$index]}"  

python sensitivity_parallel.py $blend_name
