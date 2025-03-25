#!/bin/bash

# Parse command-line arguments

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in

        --run) run="$2" ; shift ;;
        --n_seqs) n_seqs="$2" ; shift ;;
        --relax_cycles) relax_cycles="$2" ; shift ;;
        --soluble) soluble="$2" ; shift ;;
        --distance) distance="$2" ; shift ;;
        --t) t="$2" ; shift ;;
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift # Shift past the current argument
done

# Display the parsed values
source /apps/profile.d/load_all.sh
conda activate dl_binder_design

machine=`hostname`
echo "Current machine $machine"
input_silent="output/run_${run}/run_${run}_design_${t}_input.silent"
silent_out=`echo "$input_silent" | sed 's#.silent#_out.silent#'`
silent_point=`echo "$input_silent" | sed 's#.silent#_out.point#'`
fixed_residues_path="output/run_$i/fixed_residues_$t.json"
echo "Running pMPNN on $input_silent"

python3 -u /apps/rosetta/dl_binder_design/mpnn_fr/dl_interface_design_cutre.py -silent "$input_silent" -checkpoint_path "/apps/rosetta/dl_binder_design/mpnn_fr/ProteinMPNN/vanilla_model_weights/v_48_030.pt" -outsilent "$silent_out" -relax_cycles "$relax_cycles" -seqs_per_struct "$n_seqs" -checkpoint_name "$silent_point" #requires using our own modified version of mpnn, in which temp.pdb has a different name

# if [ "$soluble" = "True" ]; then
#     python3 -u /apps/scripts/protein_design/scripts/fixed_residues_sol.py --input "$input_silent" --distance "$distance"
#     python3 -u /apps/rosetta/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent "$input_silent" -checkpoint_path "/apps/rosetta/dl_binder_design/mpnn_fr/ProteinMPNN/soluble_model_weights/v_48_020.pt" -outsilent "$output_silent" -relax_cycles "$relax_cycles" -seqs_per_struct "$n_seqs" -checkpoint_name "$silent_point" -fixed_residues_path "$fixed_residues_path"
# fi
echo "done"

