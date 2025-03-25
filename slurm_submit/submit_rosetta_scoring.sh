#!/bin/bash

# Number of desired cpus:
#SBATCH --nodes=1
#SBATCH -p RFD
#SBATCH --open-mode=append
#SBATCH --exclusive=user
#SBATCH -n 48
#SBATCH -o slurm_logs/%j.out
#SBATCH -e slurm_logs/%j.err

# Activate environment
source /apps/profile.d/load_all.sh
conda activate dl_binder_design

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --folder)
            folder="$2"
            shift # Shift past the argument value
            ;;
        --run)
            run="$2"
            shift # Shift past the argument value
            ;;
    esac
    shift # Shift past the current argument
done

# Display the parsed values
machine=`hostname`
echo "Current machine $machine"
echo "Rosetta will be computed"

# Run!
python3 /apps/scripts/protein_design/scripts/rosetta_scoring_tools.py --silent $folder/run_${run}_input_out_af2.silent

echo "done"

