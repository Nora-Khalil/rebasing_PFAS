#!/bin/sh

#SBATCH --time=2-00:00:00
#SBATCH --job-name=CH3F_all_families
#SBATCH --error=rmg.slurm.log
#SBATCH --output=output.slurm.log
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --partition=short


python-jl /home/khalil.nor/Code/RMG-Py/rmg.py input.py