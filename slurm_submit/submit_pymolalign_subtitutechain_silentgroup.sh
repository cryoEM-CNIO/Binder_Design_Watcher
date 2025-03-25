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
        --pymol_template)
            pymol_template="$2"
            shift # Shift past the argument value
            ;;
        --pymol_reference)
            pymol_reference="$2"
            shift # Shift past the argument value
            ;;
        --pymol_output)
            pymol_output="$2"
            shift # Shift past the argument value
            ;;
        --sub_input)
            sub_input="$2"
            shift # Shift past the argument value
            ;;
        --sub_output)
            sub_output="$2"
            shift # Shift past the argument value
            ;;
        --silent_output)
            silent_output="$2"
            shift # Shift past the argument value
            ;;
        --rfd_ndesigns)
            rfd_ndesigns="$2"
            shift # Shift past the argument value
            ;;

        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
    shift # Shift past the current argument
done

# Display the parsed values
machine=`hostname`
echo "Current machine $machine"
echo "Pymol will align $pymol_template to $pymol_reference and save it as $pymol_output"
echo "Chain B will be substituted by the aligned $pymol_output file in all pdb files in the $sub_input and saved into $sub_output folder"
echo "4 silent files with all ready pdb files will be made"

# Run!
#run align_pymol.py
pymol -c /apps/scripts/protein_design/scripts/pymol_align.py --reference $pymol_reference --moving $pymol_template --output $pymol_output

# Substitute chains
echo "Substituting chain B for original template"
python3 /apps/scripts/protein_design/scripts/substituting_chains.py --input-dir $sub_input --chain-to-substitute B --chain-to-remove B --second-structure $pymol_output --output-dir $sub_output 

#run silent to split in 4 groups each for each gpu
echo "Creating silent files"
cd $sub_output ; /apps/rosetta/dl_binder_design/include/silent_tools/silentfrompdbsparallel *substituted.pdb > "$silent_output"
# echo "Splitting pdbs into 4 split files"
# /apps/rosetta/dl_binder_design/include/silent_tools/silentsplitshuf "$silent_output" $((rfd_ndesigns / 4))

echo "done"

