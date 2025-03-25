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

: '
Submit the pMPNN soluble for the residues that are not involved in the interaction
in order to improve the solubility of the designs

Since microruns works like a clock and name changes is going to mess all up, I am gonna overwritte the silent although that should not be done
'

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --input_silent)
            input_pmp="$2"
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
        --fixed_residues_path)
            fixed_residues="$2"
            shift
            ;;
        --distance)
            distance="$2"
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
echo "Running pMPNN on $input_pmp"

silent_point=`echo "$input_pmp" | sed s'#.silent#_sol.point#'`
input_silent=`echo "$input_pmp" | sed s',.silent,_out.silent,'`
output_silent=`echo "$input_pmp" | sed s',.silent,_sol_out.silent,'`

echo "$input_silent"

python3 -u /apps/scripts/protein_design/scripts/fixed_residues_sol.py --input "$input_silent" --distance "$distance"
python3 -u /apps/rosetta/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent "$input_silent" -checkpoint_path "/apps/rosetta/dl_binder_design/mpnn_fr/ProteinMPNN/soluble_model_weights/v_48_020.pt" -outsilent "$output_silent" -relax_cycles "$relax_cycles" -seqs_per_struct "$n_seqs" -checkpoint_name "$silent_point" -fixed_residues_path "$fixed_residues"

echo "done"

