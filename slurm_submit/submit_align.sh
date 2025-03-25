#!/bin/bash

# Number of desired cpus:
#SBATCH --nodes=1
#SBATCH -p RFD
#SBATCH --open-mode=append
#SBATCH --exclusive=user
#SBATCH -o slurm_logs/%j.out
#SBATCH -e slurm_logs/%j.err

# Activate environment
source /apps/profile.d/load_all.sh
conda activate dl_binder_design

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --template) template="$2" ; shift ;;
        --run) run="$2" ; shift ;;
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift 
done

# Variables declaring
input_dir="output/run_${run}"
silent_output="${input_dir}/run_${run}_input.silent"
reference="output/run_${run}/run_${run}_design_0.pdb"


# Display the parsed values
machine=`hostname`
echo "Current machine $machine"
echo "Template aligning and chain substitution using Biopython for folder $input_dir"
echo "Silent output for next step will be saved as $silent_output at $input_dir"

# Run!
#run align
python3 /apps/scripts/protein_design/scripts/biopython_align.py --reference $reference --template $template --chain "B" --input_dir $input_dir  
#run silent to split in 4 groups each for each gpu
echo "Creating silent files"
/apps/rosetta/dl_binder_design/include/silent_tools/silentfrompdbsparallel "${input_dir}/*substituted.pdb" > "$silent_output"

echo "done"

