
#!/bin/bash --login
#SBATCH --job-name=mass_fraction_restart
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=71:59:00
#SBATCH --output=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Outputs/%A_%a-restart.out
#SBATCH --error=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Errors/%A_%a-restart.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monkhayd@msu.edu

set -euo pipefail

ROOT="/scratch/group/p.phy260085.000/Mass_Fraction"
JOB_ID="2215154"

PYTHON_SCRIPT="$ROOT/Inputs/N_Body_Scripts/restart_mass_fraction_simple.py"

module purge
module load Anaconda3/2022.10
source activate comets

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Original job ID: $JOB_ID"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "MPI ranks: $SLURM_NTASKS"
echo "Start: $(date)"

mpirun -np "$SLURM_NTASKS" \
    "$CONDA_PREFIX/bin/python" -u \
    "$PYTHON_SCRIPT" "$JOB_ID" "$SLURM_ARRAY_TASK_ID"

echo "Finished: $(date)"
