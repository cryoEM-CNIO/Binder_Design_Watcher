#!/bin/bash

# Number of desired cpus:
#SBATCH --nodes=1
#SBATCH -p RFD
#SBATCH --open-mode=append
#SBATCH --exclusive=user
#SBATCH --gres=gpu:1
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
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
    shift # Shift past the current argument
done

source /apps/profile.d/load_all.sh
conda activate dl_binder_design

machine=`hostname`
echo "Current machine $machine"


af2_in=`echo "$input_silent" | sed 's#\.silent#_out.silent#'`
af2_out=`echo "$input_silent" | sed 's#\.silent#_out_af2.silent#'`
af2_score=`echo "$input_silent" | sed 's#\.silent#_out_af2.sc#'`
af2_point=`echo "$input_silent" | sed 's#\.silent#_out_af2.point#'`
af2_json=$(echo "$input_silent" | sed -E 's#run_[0-9]+_input.*\.silent#pae#')
echo "Running AF2 on $af2_in"

python3 -u /apps/rosetta/dl_binder_design/af2_initial_guess/predict_cutre.py -silent "$af2_in" -outsilent "$af2_out" -scorefilename "$af2_score"  -checkpoint_name "$af2_point" -jsonfilename "$af2_json"

echo "done"
