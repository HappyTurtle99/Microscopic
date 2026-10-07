#!/bin/bash
#SBATCH --job-name=micro
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --partition=nariman
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs

export JULIA_SCRATCH_DIR="$HOME/.julia/scratch"
export JULIA_DEPOT_PATH="$HOME/.julia"

# Run the script using the Microscopic project directory
julia --project=/home/beq/gtucci/Work/DiscreteChemo/DiscreteChemo/Microscopic /home/beq/gtucci/Work/DiscreteChemo/DiscreteChemo/Microscopic/rough_shape18jun_nariman.jl 100 5.0 0.000001