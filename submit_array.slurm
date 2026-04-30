#!/bin/bash
#SBATCH --job-name=micro_array
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --partition=nariman
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-3

set -e

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs

JULIA=/home/beq/gtucci/.juliaup/bin/julia

# Parameter lists
Nsites_vals=(30.0 40.0 50.0 60.0)
Dc_vals=(0.2 1 2 6)
lambda_vals=(1.0)

# Get SLURM index
i=$SLURM_ARRAY_TASK_ID

# Map index → parameters
n_idx=$(( i % ${#Nsites_vals[@]} ))
d_idx=$(( (i / ${#Nsites_vals[@]}) % ${#Dc_vals[@]} ))
l_idx=$(( i / (${#Nsites_vals[@]} * ${#Dc_vals[@]}) ))

Nsites=${Nsites_vals[$n_idx]}
Dc=${Dc_vals[$d_idx]}
lambda1=${lambda_vals[$l_idx]}

echo "Running job $i with:"
echo "Nsites=$Nsites Dc=$Dc lambda1=$lambda1"

# Run simulation (adjust argument order if needed)
$JULIA --project=. simulation_1sp.jl $Nsites $Dc $lambda1