#!/bin/bash --login

#SBATCH --job-name=mass_fraction_sim          # Name of Job. Its whatever you want to call it.  
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1          
#SBATCH --ntasks=46                      #Dang suggest 128?    # Number of tasks. This is also the number of cores  # SLURM defaults to 1 but we specify anyway

## Use nodes keyword
## also ntasks-per-node
#SBATCH --mem=1G                              # Memory per node   # Specify "M" or "G" for MB and GB respectively
#SBATCH --time=9:59:00                         # Wall time         # Format: "minutes", "hours:minutes:seconds",      # "days-hours", or "days-hours:minutes"
#SBATCH --output=/scratch/group/p.phy260085.000/mass_fraction_hayden/Outputs/SLURM_Outputs/%A-SLURM_output.out                # This will save your outputs to a file  %x: job name, %j: job ID
#SBATCH --error=/scratch/group/p.phy260085.000/mass_fraction_hayden/Outputs/SLURM_Outputs/%A-SLURM_error.out

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


conda info --envs
srun which python
#srun bash -c 'python -u /scratch/group/p.phy260085.000/mass_fraction_hayden/Inputs/N_Body_Scripts/local_disc_aces.py "$1" "$2" "$BASHPID" "$3"' _ \
 #   "$SLURM_ARRAY_JOB_ID" \
  #  "$SLURM_ARRAY_TASK_ID" \
   # "$PARAM_FILE"

srun bash -c "\"$PYTHON_BIN\" -u /scratch/group/p.phy260085.000/mass_fraction_hayden/Inputs/N_Body_Scripts/aces_local_disc_slim_output.py \"\$1\" \"\$2\" \"\$BASHPID\" \"\$3\"" _ \
    "$SLURM_ARRAY_JOB_ID" \
    "$SLURM_ARRAY_TASK_ID" \
    "$PARAM_FILE"


