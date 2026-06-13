#!/bin/bash
#SBATCH --job-name=Dckappasweep
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --partition=rostam
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-3

set -e

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs

JULIA=/home/beq/gtucci/.juliaup/bin/julia

# As of Jun 3rd 18th, Nsites means Dc and Dc means kappa
Nsites_vals=(0.01 0.075 0.3 1.2)
Dc_vals=(0.025)

# Get SLURM index
i=$SLURM_ARRAY_TASK_ID

# Map index → parameters
n_idx=$(( i % ${#Nsites_vals[@]} ))
d_idx=$(( (i / ${#Nsites_vals[@]}) % ${#Dc_vals[@]} ))

Nsites=${Nsites_vals[$n_idx]}
Dc=${Dc_vals[$d_idx]}

echo "Running job $i with:"
echo "Dc=$Nsites kappa=$Dc"

$JULIA Dc_sweep_13Jun_Dc.jl $Nsites $Dc
