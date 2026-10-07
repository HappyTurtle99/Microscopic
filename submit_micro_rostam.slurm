#!/bin/bash
#SBATCH --job-name=micro
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --partition=rostam
#SBATCH --time=165:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs

export JULIA_SCRATCH_DIR="$HOME/.julia/scratch"
export JULIA_DEPOT_PATH="$HOME/.julia"

#REGULAR
# julia --project=/home/beq/gtucci/Work/DiscreteChemo/DiscreteChemo/Microscopic/home/beq/gtucci/Work/DiscreteChemo/DiscreteChemo/Microscopic/rough_shape18jun_rostam.jl 170 1.0 150.0

#CONTINUE
julia --project=/home/beq/gtucci/Work/DiscreteChemo/DiscreteChemo/Microscopic /home/beq/gtucci/Work/DiscreteChemo/DiscreteChemo/Microscopic/continue_sim.jl