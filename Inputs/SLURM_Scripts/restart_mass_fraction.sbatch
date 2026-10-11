#!/bin/bash --login
#SBATCH --job-name=mass_fraction_restart
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=71:59:00
#SBATCH --error=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Errors/%A_%a-restart_error.out
#SBATCH --output=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Outputs/%A_%a-restart_output.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monkhayd@msu.edu

set -euo pipefail

ROOT="/scratch/group/p.phy260085.000/Mass_Fraction"
PYTHON_SCRIPT="$ROOT/Inputs/N_Body_Scripts/restart_mass_fraction.py"
SLURM_SCRIPT="$ROOT/Inputs/SLURM_Scripts/restart_mass_fraction.sbatch"
ARCHIVE_ROOT="$ROOT/Outputs/Sim_Archives"
JOB_ID="2199100"
MAX_CONCURRENT_MASSES=3

# When invoked as `bash restart_mass_fraction.sbatch` or
# `sbatch restart_mass_fraction.sbatch` without --array, discover the
# mass directories and submit this same file as an array job.
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    if [[ -z "$JOB_ID" || "$JOB_ID" == "REPLACE_WITH_ORIGINAL_JOB_ID" ]]; then
        echo "Set JOB_ID near the top of this SLURM script." >&2
        exit 1
    fi

    shopt -s nullglob
    mass_dirs=("$ARCHIVE_ROOT/${JOB_ID}_mass_"*)
    count=0
    for dir in "${mass_dirs[@]}"; do
        if [[ -d "$dir" && "${dir##*_mass_}" =~ ^[0-9.eE+-]+$ ]]; then
            count=$((count + 1))
        fi
    done
    if ((count == 0)); then
        echo "No mass archive directories found for job $JOB_ID" >&2
        exit 1
    fi

    mkdir -p "$ROOT/Outputs/SLURM_Outputs/Errors" \
             "$ROOT/Outputs/SLURM_Outputs/Outputs"
    echo "Submitting $count mass tasks ($MAX_CONCURRENT_MASSES concurrent), original job: $JOB_ID"
    sbatch --array="0-$((count-1))%${MAX_CONCURRENT_MASSES}" "$SLURM_SCRIPT"
    exit 0
fi

# Array mode: exactly one mass, selected by its sorted directory index.
export SLURM_EXPORT_ENV=ALL
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module purge
module load Anaconda3/2022.10
source activate comets
PYTHON_BIN="$CONDA_PREFIX/bin/python"

echo "=========================================="
echo "Mass Fraction Restart"
echo "Date:             $(date)"
echo "Host:             $(hostname)"
echo "Array job ID:     ${SLURM_ARRAY_JOB_ID:-unknown}"
echo "Original JOB_ID:  $JOB_ID"
echo "Mass index:       $SLURM_ARRAY_TASK_ID"
echo "MPI ranks:        $SLURM_NTASKS"
echo "Python:           $PYTHON_BIN"
echo "=========================================="

mpirun -np "$SLURM_NTASKS" \
    "$PYTHON_BIN" -u "$PYTHON_SCRIPT" "$JOB_ID" "$SLURM_ARRAY_TASK_ID"
