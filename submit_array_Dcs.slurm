#!/bin/bash
#SBATCH --job-name=santi_dc
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --partition=rostam
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-4

set -e

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs

JULIA=/home/beq/gtucci/.juliaup/bin/julia

# As of May 18th, Nsites means Nsites and Dc means Tfinal
Nsites_vals=(10.0 15.0 20.0 25.0 45.0)
Dc_vals=(1)

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
$JULIA continue_sim.jl $Nsites $Dc
