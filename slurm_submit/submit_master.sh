#!/bin/bash

# Number of desired cpus:
#SBATCH --nodes=1
#SBATCH -p RFD
#SBATCH --open-mode=append
#SBATCH --gres=gpu:1
#SBATCH --exclusive
#SBATCH --cpus-per-gpu=12
#SBATCH -o slurm_logs/%j.out
#SBATCH -e slurm_logs/%j.err

## Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input) input="$2" ; shift ;;
        -t|--template) template="$2" ; shift ;;
        -r|--run) run="$2" ; shift ;;
        -c|--rfd_contigs) rfd_contigs="$2" ; shift  ;;    
        -h|--rfd_hotspots) rfd_hotspots="$2" ; shift  ;;    
        -nd|--rfd_ndesigns) rfd_ndesigns="$2" ; shift  ;;    
        -np|--pmp_nseqs) pmp_nseqs="$2" ; shift  ;;    
        -rc|--pmp_relax_cycles) pmp_relax_cycles="$2" ; shift  ;;   
        -pd|--partial_diff) partial_diff="$2" ; shift  ;; 
        -nst|--noise_steps) noise_steps="$2" ; shift  ;;
        -nsc|--noise_scale) noise_scale="$2" ; shift  ;;
        -ck|--ckp) ckp="$2" ; shift  ;; #Add the path to the checkpoint to add weight toward some fold
        -w|--node) node="$2" ; shift  ;; # Provide a specific name of a node to submit to this node with -w. If not provided it will be as usual.
        -sp|--soluble_pMPNN) soluble="$2" ; shift  ;; #Add a pMPNNsol step to the non interacting interface to increase the solubility
        -d|--distance) distance="$2" ; shift  ;; # Distance to fix residues in the interface
        -hn|--hits_number) hits_number="$2" ; shift ;;
        -cr|--core) core="$2" ; shift ;; # proportion of core residues to make the filtering
        -re|--residues) residues="$2" ; shift ;; # Residues to fix, useful for scaffolding
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift  # Shift past the current argument
done

# Get available GPUs from SLURM
GPUS_AVAILABLE=$(nvidia-smi --query-gpu=index --format=csv,noheader | tr '\n' ' ')
echo "GPUs available: $GPUS_AVAILABLE"

t=1

for GPU_ID in $GPUS_AVAILABLE; do
    echo "Using $GPU_ID"
    (
        export CUDA_VISIBLE_DEVICES=$GPU_ID
        LOG_DIR="slurm_logs/${SLURM_JOB_ID}_gpu${GPU_ID}"
        mkdir -p "$LOG_DIR"
        # ----------------------------------------
        # 1. RFD
        # ----------------------------------------

        bash /apps/scripts/protein_design/master_scripts/rfd.sh \
            --run "$run" --t "$t" \
            --input_pdb "$input" --contigmap_descriptor "$rfd_contigs" \
            --designs_n "$rfd_ndesigns" --noise_steps "$noise_steps" \
            --noise_scale "$noise_scale" --ckp "$ckp" \
            --partial_diff "$partial_diff" --residues "$residues" --hotspots_descriptor "$rfd_hotspots"  > "$LOG_DIR/rfd.out" 2> "$LOG_DIR/rfd.err"
        wait

        # --------------------------------------------
        # 2. Aligning + filtering
        # --------------------------------------------

        bash /apps/scripts/protein_design/master_scripts/aligning_filtering.sh --template "$template" --run "$run" --t "$t" --core "$core" --residues "$residues" > "$LOG_DIR/aligning_filtering.out" 2> "$LOG_DIR/aligning_filtering.err"
        wait
        # --------------------------------------------
        # 3 pMPNN
        # --------------------------------------------

        bash /apps/scripts/protein_design/master_scripts/pmpnn.sh --run "$run" --t "$t" --n_seqs "$pmp_nseqs" --relax_cycles "$pmp_relax_cycles" --soluble "$soluble" --distance "$distance" > "$LOG_DIR/pmpnn.out" 2> "$LOG_DIR/pmpnn.err"
        wait
        
        # --------------------------------------------
        # 4 Scoring(AF2IG + PyRosetta)
        # --------------------------------------------

        bash /apps/scripts/protein_design/master_scripts/scoring.sh --run "$run" --t "$t" > "$LOG_DIR/scoring.out" 2> "$LOG_DIR/scoring.err"
        wait
    ) &
    ((t=t+1))
done
wait
# --------------------------------------------
# 5 Finish Microrun
# --------------------------------------------

bash /apps/scripts/protein_design/master_scripts/ending.sh --number "$hits_number" --run "$run"