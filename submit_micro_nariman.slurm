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

echo "copy 2: lambda 5 Tfinal 75"

/home/beq/gtucci/.juliaup/bin/julia 
# continue_sim.jl
rough_shape18jun_hom.jl 200 5.0 0.001