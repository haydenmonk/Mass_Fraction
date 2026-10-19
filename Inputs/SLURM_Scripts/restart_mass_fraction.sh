#!/bin/bash --login
#SBATCH --job-name=mass_fraction_restart
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=71:59:00
#SBATCH --output=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Outputs/%j-restart.out
#SBATCH --error=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Errors/%j-restart.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monkhayd@msu.edu

set -euo pipefail

ROOT="/scratch/group/p.phy260085.000/Mass_Fraction"
JOB_ID="2215154"

ARCHIVE_ROOT="$ROOT/Outputs/Sim_Archives"
PYTHON_SCRIPT="$ROOT/Inputs/N_Body_Scripts/restart_mass_fraction.py"

module purge
module load Anaconda3/2022.10
source activate comets

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Expected directories:
# Outputs/Sim_Archives/2215154_mass_0.0001/
# Outputs/Sim_Archives/2215154_mass_0.001/
# Outputs/Sim_Archives/2215154_mass_0.01/
# ... (actual values are discovered automatically)

shopt -s nullglob
MASS_DIRS=("$ARCHIVE_ROOT/${JOB_ID}_mass_"*)

echo "Found ${#MASS_DIRS[@]} candidate mass directories."

if (( ${#MASS_DIRS[@]} == 0 )); then
    echo "ERROR: No mass directories found."
    exit 1
fi

# Restart one mass at a time, using 24 MPI ranks per mass.
for i in "${!MASS_DIRS[@]}"; do
    echo "=================================="
    echo "Mass index: $i"
    echo "Directory:  ${MASS_DIRS[$i]}"
    echo "Date:       $(date)"
    echo "=================================="

    mpirun -np "$SLURM_NTASKS" \
        "$CONDA_PREFIX/bin/python" -u \
        "$PYTHON_SCRIPT" "$JOB_ID" "$i"
done

echo "All mass restarts completed."

