#!/bin/sh

#SBATCH --time=20-00:00:00
#SBATCH --job-name=S_CH3F
#SBATCH --error=rmg.slurm.log
#SBATCH --output=output.slurm.log
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --partition=west


python-jl /home/khalil.nor/Code/RMG-Py/rmg.py input.py