#!/bin/bash
#SBATCH --job-name=SShydro
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --partition=rostam
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-8

set -e

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs

JULIA=/home/beq/gtucci/.juliaup/bin/julia

# As of Jun 25th, Nsites means 1/ h and Dc is a decoy
Nsites_vals=(75 100 150 175 200 225 250 300 325)
Dc_vals=(1.0)

# Get SLURM index
i=$SLURM_ARRAY_TASK_ID

# Map index → parameters
n_idx=$(( i % ${#Nsites_vals[@]} ))
d_idx=$(( (i / ${#Nsites_vals[@]}) % ${#Dc_vals[@]} ))

Nsites=${Nsites_vals[$n_idx]}
Dc=${Dc_vals[$d_idx]}

echo "Running job $i with:"
echo "1 / h =$Nsites kappa=$Dc"

$JULIA rough_shape18jun.jl $Nsites $Dc
