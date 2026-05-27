#!/bin/bash
#SBATCH --job-name=micro_array
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --partition=rostam
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-20

set -e

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs

JULIA=/home/beq/gtucci/.juliaup/bin/julia

# As of May 18th, Nsites means Nsites and Dc means Tfinal
Nsites_vals=(1500 2200 2900 3600 4200 4900)
Dc_vals=(50.0 100.0 150.0)

# Get SLURM index
i=$SLURM_ARRAY_TASK_ID

# Map index → parameters
n_idx=$(( i % ${#Nsites_vals[@]} ))
d_idx=$(( (i / ${#Nsites_vals[@]}) % ${#Dc_vals[@]} ))

Nsites=${Nsites_vals[$n_idx]}
Dc=${Dc_vals[$d_idx]}

echo "Running job $i with:"
echo "Nsites=$Nsites T=$Dc"

# Run simulation (adjust argument order if needed)
$JULIA window_tests18May2026.jl $Nsites $Dc