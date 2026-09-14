#!/bin/bash --login

#SBATCH --job-name=mass_fraction_mpi

# Six active masses: indices 0-5
# %1 means only one mass runs at a time
#SBATCH --array=0-5%1

# Number of MPI ranks working on EACH mass
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1

#SBATCH --mem=10G
#SBATCH --time=12:59:00

#SBATCH --error=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Errors/%A_%a-SLURM_error.out
#SBATCH --output=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Outputs/%A_%a-SLURM_output.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=monkhayd@msu.edu


set -euo pipefail


# ============================================================
# Configuration
# ============================================================

PYTHON_SCRIPT="/scratch/group/p.phy260085.000/Mass_Fraction/Inputs/N_Body_Scripts/aces_local_disc_mpi_multi.py"



PARAM_FILE="unused"


# ============================================================
# Environment
# ============================================================

export SLURM_EXPORT_ENV=ALL

# Prevent NumPy / BLAS from spawning additional threads
# inside every MPI rank.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module purge
module load Anaconda3/2022.10

source activate comets

PYTHON_BIN="$CONDA_PREFIX/bin/python"


# ============================================================
# Diagnostics
# ============================================================

echo "=========================================="
echo "Mass Fraction MPI Run"
echo "=========================================="
echo "Date:               $(date)"
echo "Host:               $(hostname)"
echo "Array job ID:       $SLURM_ARRAY_JOB_ID"
echo "Array task ID:      $SLURM_ARRAY_TASK_ID"
echo "MPI ranks:          $SLURM_NTASKS"
echo "CPUs per task:      $SLURM_CPUS_PER_TASK"
echo "Nodes allocated:    $SLURM_JOB_NUM_NODES"
echo "Python:             $PYTHON_BIN"
echo "Conda environment:  $CONDA_PREFIX"
echo "=========================================="
echo


# ============================================================
# Launch MPI Python program
# ============================================================

srun --kill-on-bad-exit=1 \
    "$PYTHON_BIN" -u "$PYTHON_SCRIPT" \
    "$SLURM_ARRAY_JOB_ID" \
    "$SLURM_ARRAY_TASK_ID" \
    "$PARAM_FILE"
