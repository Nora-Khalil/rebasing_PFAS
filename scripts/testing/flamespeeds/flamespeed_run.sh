#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --job-name=fcS_CH3F
#SBATCH --error=fc.slurm.log
#SBATCH --output=fc_output.slurm.log
#SBATCH --partition=short 

python flamespeeds.py '/work/westgroup/nora/Code/projects/halogens/refrigerants/halogens_paper/changing_RMG_Py/on_USNCM_blends/singles/CH3F/cantera/chem_annotated.yaml' 'CH3F(1)'
