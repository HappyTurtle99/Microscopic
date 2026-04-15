#!/bin/bash
#SBATCH --job-name=micro
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --partition=nariman
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

cd $SLURM_SUBMIT_DIR

source /home/beq/gtucci/miniconda3/etc/profile.d/conda.sh
conda activate base

python simulation_1sp.py 25 0.8 100
