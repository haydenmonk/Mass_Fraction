#!/bin/bash --login

#SBATCH --job-name=mpi_test
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --output=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Outputs/mpi_test_%j.out
#SBATCH --error=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Outputs/mpi_test_%j.err

set -euo pipefail

module purge
module load Anaconda3/2022.10

source activate comets

PYTHON_BIN="$CONDA_PREFIX/bin/python"

echo "Python: $PYTHON_BIN"
echo "SLURM_NTASKS: $SLURM_NTASKS"

mpirun -np "$SLURM_NTASKS" \
    "$PYTHON_BIN" -u /scratch/group/p.phy260085.000/Mass_Fraction/Inputs/N_Body_Scripts/mpi_test.py
