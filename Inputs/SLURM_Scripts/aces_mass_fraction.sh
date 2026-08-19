#!/bin/bash --login

#SBATCH --job-name=mass_fraction_sim          # Name of Job. Its whatever you want to call it.  
#SBATCH --array=0-80
#SBATCH --cpus-per-task=1          
#SBATCH --ntasks=500                      #Dang suggest 128?    # Number of tasks. This is also the number of cores  # SLURM defaults to 1 but we specify anyway

## Use nodes keyword
## also ntasks-per-node
#SBATCH --mem=100G                              # Memory per node   # Specify "M" or "G" for MB and GB respectively
#SBATCH --time=5:59:00                         # Wall time         # Format: "minutes", "hours:minutes:seconds",      # "days-hours", or "days-hours:minutes"
#Sbatch --error=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Errors/%A-SLURM_error.out
#SBATCH --output=/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/SLURM_Outputs/Outputs/%A-SLURM_output.out
# There next two are optional and have been commented out. Uncomment if you want emails to send to you about your job 
#SBATCH --mail-type=ALL                       # Mail type         # e.g., which events trigger email notifications
#SBATCH --mail-user=monkhayd@msu.edu        # Mail address 

# If you load up a conda environment this is how to do that
export SLURM_EXPORT_ENV=ALL
module purge 
module load Anaconda3/2022.10
source activate comets
PYTHON_BIN="$CONDA_PREFIX/bin/python"
echo "Using interpreter: $PYTHON_BIN"
which python
echo $CONDA_PREFIX
ls $CONDA_PREFIX/bin | grep python
echo "Starting mass fraction simulations"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Nodes allocated: $SLURM_JOB_NUM_NODES"
echo "CPUs allocated on this node: $SLURM_CPUS_ON_NODE"
echo "Working directory: $(pwd)"
echo

#srun bash -c 'python -u /scratch/group/p.phy260085.000/mass_fraction_hayden/Inputs/N_Body_Scripts/local_disc_aces.py "$1" "$2" "$BASHPID" "$3"' _ \
 #   "$SLURM_ARRAY_JOB_ID" \
  #  "$SLURM_ARRAY_TASK_ID" \
   # "$PARAM_FILE"

srun bash -c "\"$PYTHON_BIN\" -u /scratch/group/p.phy260085.000/Mass_Fraction/Inputs/N_Body_Scripts/aces_local_disc_analytic.py \"\$1\" \"\$2\" \"\$BASHPID\" \"\$3\"" _ \
    "$SLURM_ARRAY_JOB_ID" \
    "$SLURM_ARRAY_TASK_ID" \
    "$PARAM_FILE"


