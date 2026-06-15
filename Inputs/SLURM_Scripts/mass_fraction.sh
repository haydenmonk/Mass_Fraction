#!/bin/bash --login

#SBATCH --job-name=mass_fraction_sim          # Name of Job. Its whatever you want to call it.  
#SBATCH --array=0-160
#SBATCH --cpus-per-task=1          
#SBATCH --ntasks=200                      #Dang suggest 128?    # Number of tasks. This is also the number of cores  # SLURM defaults to 1 but we specify anyway

## Use nodes keyword
## also ntasks-per-node
#SBATCH --mem=10G                              # Memory per node   # Specify "M" or "G" for MB and GB respectively
#SBATCH --time=47:59:00                         # Wall time         # Format: "minutes", "hours:minutes:seconds",      # "days-hours", or "days-hours:minutes"
#SBATCH --output=/mnt/home/monkhayd/Simulations/Ejection_modeling/Mass_Fraction/Outputs/SLURM_Outputs/Outputs/%A-SLURM_output.out                # This will save your outputs to a file  %x: job name, %j: job ID
#SBATCH --error=/mnt/home/monkhayd/Simulations/Ejection_modeling/Mass_Fraction/Outputs/SLURM_Outputs/Errors/%A-SLURM_error.out

# There next two are optional and have been commented out. Uncomment if you want emails to send to you about your job 
#SBATCH --mail-type=ALL                       # Mail type         # e.g., which events trigger email notifications
#SBATCH --mail-user=monkhayd@msu.edu        # Mail address 

# If you load up a conda environment this is how to do that
module purge 
. ~/Software/miniconda3/etc/profile.d/conda.sh

conda activate comets


PARAM_FILE="/mnt/home/monkhayd/Simulations/Ejection_modeling/Mass_Fraction/Inputs/Parameter Files/Parameter_CSVs/m_a_parameter_values.csv"

srun bash -c 'python -u /mnt/home/monkhayd/Simulations/Ejection_modeling/Mass_Fraction/Inputs/N_Body_Scripts/local_disc.py "$1" "$2" "$BASHPID" "$3"' _ \
    "$SLURM_ARRAY_JOB_ID" \
    "$SLURM_ARRAY_TASK_ID" \
    "$PARAM_FILE"