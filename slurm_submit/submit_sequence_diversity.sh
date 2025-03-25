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

# Display the parsed values
source /apps/profile.d/load_all.sh
conda activate dl_binder_design

hits_number=10000
## Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -r|--run) run="$2" ; shift ;;
        -ns|--nseqs) pmp_nseqs="$2" ; shift  ;;    
        -fr|--fr) pmp_relax_cycles="$2" ; shift  ;;   
        -hn|--hits_number) hits_number="$2" ; shift ;;
        -f|--fixed) fixed="$2" ; shift ;; 
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift  # Shift past the current argument
done

# Get available GPUs from SLURM
GPUS_AVAILABLE=$(nvidia-smi --query-gpu=index --format=csv,noheader | tr '\n' ' ')
echo "GPUs available: $GPUS_AVAILABLE"

t=1

for GPU_ID in $GPUS_AVAILABLE; do
    echo $GPU_ID
    (
        export CUDA_VISIBLE_DEVICES=$GPU_ID
        LOG_DIR="slurm_logs/${SLURM_JOB_ID}_gpu${GPU_ID}"
        mkdir -p "$LOG_DIR"

        # --------------------------------------------
        # 1 Generate the silent file
        # --------------------------------------------

        /apps/rosetta/dl_binder_design/include/silent_tools/silentrename initial_input.silent "run_${run}_design_${t}" > "output/run_${run}/run_${run}_design_${t}_input.silent" 
        wait
        # --------------------------------------------
        # 2 pMPNN
        # --------------------------------------------

        bash /apps/scripts/protein_design/master_scripts/pmpnn.sh --run "$run" --t "$t" --n_seqs "$pmp_nseqs" --relax_cycles "$pmp_relax_cycles" --soluble "$soluble" --distance "$distance" > "$LOG_DIR/pmpnn.out" 2> "$LOG_DIR/pmpnn.err"
        wait
        # --------------------------------------------
        # 3 Scoring(AF2IG + PyRosetta)
        # --------------------------------------------

        bash /apps/scripts/protein_design/master_scripts/scoring.sh --run "$run" --t "$t" > "$LOG_DIR/scoring.out" 2> "$LOG_DIR/scoring.err"
        wait
    ) &
    ((t=t+1))
done
wait
# --------------------------------------------
# 4 Finish Microrun
# --------------------------------------------

bash /apps/scripts/protein_design/master_scripts/ending.sh --number "$hits_number" --run "$run"