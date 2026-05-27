#!/bin/bash
#SBATCH --job-name=micro
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --partition=rostam
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs

/home/beq/gtucci/.juliaup/bin/julia window_tests18May2026.jl 4900 0.001
