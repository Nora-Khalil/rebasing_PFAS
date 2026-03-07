#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=6-00:00:00
#SBATCH --job-name=autoTrue
#SBATCH --error=rmg.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --partition=west

python /home/khalil.nor/Code/RMG-Py/rmg.py input.py

