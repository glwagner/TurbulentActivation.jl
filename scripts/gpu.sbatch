#!/bin/bash
#SBATCH -J turbulent-activation
#SBATCH -p gpua100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=6
#SBATCH --gpus=1
#SBATCH --mem=0
#SBATCH -t 12:00:00
#SBATCH -o /home/greg_aeolus_earth/TurbulentActivation.jl/runs/slurm-%j.out
# Usage: sbatch scripts/gpu.sbatch scripts/run_droplets.jl 64 64 32 5 10000 parity
mkdir -p $HOME/TurbulentActivation.jl/runs && cd $HOME/TurbulentActivation.jl/runs
script=$1; shift
julia --startup-file=no $JULIA_FLAGS --project=$HOME/TurbulentActivation.jl $HOME/TurbulentActivation.jl/$script "$@"
