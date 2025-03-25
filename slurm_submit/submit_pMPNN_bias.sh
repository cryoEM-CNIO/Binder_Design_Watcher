#!/bin/bash

# Number of desired cpus:
#SBATCH --nodes=1
#SBATCH -p RFD
#SBATCH --open-mode=append
#SBATCH --gres=gpu:1
#SBATCH --exclusive=user
#SBATCH --cpus-per-gpu=12
#SBATCH -o slurm_logs/%j.out
#SBATCH -e slurm_logs/%j.err

# Parse command-line arguments

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --input_silent)
            input_silent="$2"
            shift # Shift past the argument value
            ;;
        --n_seqs)
            n_seqs="$2"
            shift
            ;;
        --relax_cycles)
            relax_cycles="$2"
            shift
            ;;
        --bias)
            bias="$2"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
    shift # Shift past the current argument
done

# Display the parsed values
source /apps/profile.d/load_all.sh
conda activate dl_binder_design

machine=`hostname`
echo "Current machine $machine"
echo "Running pMPNN on $input_silent"

silent_out=`echo "$input_silent" | sed 's#.silent#_out.silent#'`
silent_point=`echo "$input_silent" | sed 's#.silent#_out.point#'`

python3 -u /apps/rosetta/dl_binder_design/mpnn_fr/dl_interface_design.py -silent "$input_silent" -checkpoint_path "/apps/rosetta/dl_binder_design/mpnn_fr/ProteinMPNN/vanilla_model_weights/v_48_030.pt" -outsilent "$silent_out" -relax_cycles "$relax_cycles" -seqs_per_struct "$n_seqs" -checkpoint_name "$silent_point" -bias_AA_jsonl "$bias"

echo "done"

