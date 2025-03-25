#!/bin/bash

# Number of desired cpus:
#SBATCH --nodes=1
#SBATCH -p RFD
#SBATCH --open-mode=append
#SBATCH --exclusive=user
#SBATCH -o slurm_logs/%j.out
#SBATCH -e slurm_logs/%j.err

# Prepare
source /apps/profile.d/load_all.sh

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --number) number="$2" ; shift ;;
        --output) output="$2" ; shift ;;
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift # Shift past the current argument
done

# Display the parsed values
machine=`hostname`
echo "Current machine $machine"

echo "input_silent: $output"

# Run
touch $output
python3 /apps/scripts/protein_design/scripts/run_ending.py --number "$number"

echo "Touched $output"
