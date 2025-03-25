#!/bin/bash

# Number of desired cpus:
#SBATCH --nodes=1
#SBATCH -p RFD
#SBATCH --open-mode=append
#SBATCH --gres=gpu:2
#SBATCH --exclusive
#SBATCH --cpus-per-gpu=12
#SBATCH -o %j.out
#SBATCH -e %j.err

# Function to display usage information
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h, --help                Display this help message"
    echo "  -i, --input <x.pdb>       Input structure for partial diffusion"
    echo "  -n, --number N            Produce N designs"
    echo "  -nst, --noise_steps N     Noise steps"
    echo "  -nsc, --noise_scale N     Noise scale"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            ;;
        -i|--input)
            if [[ $# -lt 2 ]] || [[ $2 == -* ]]; then
                echo "Error: --input requires an argument"
                usage
            fi
            input_pdb="$2"
            shift # past argument
            shift # past value
            ;;
        -n|--number)
            if [[ $# -lt 2 ]] || [[ $2 == -* ]]; then
                echo "Error: --number requires an argument"
                usage
            fi
            designs_n="$2"
            shift # past argument
            shift # past value
            ;;
        -nst|--noise_steps)
            if [[ $# -lt 2 ]] || [[ $2 == -* ]]; then
                echo "Error: --noise_steps requires an argument"
                usage
            fi
            noise_steps="$2"
            shift # past argument
            shift # past value
            ;;
        -nsc|--noise_scale)
            if [[ $# -lt 2 ]] || [[ $2 == -* ]]; then
                echo "Error: --noise_scale requires an argument"
                usage
            fi
            noise_scale="$2"
            shift # past argument
            shift # past value
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Activate environment
source /apps/profile.d/load_all.sh
conda activate SE3nv4090

# Display parsed values
echo "Parsed command line arguments:"
echo "Input PDB for partially diffuse: $input_pdb"
echo "Number of designs: $designs_n"
echo "Noising steps: $noise_steps"
echo "Noising scale: $noise_scale"

# Get contigs map from PDB
rfd_contigs=$(contigs_map_getter_pd.py --input $input_pdb)
echo "Contigs map: $rfd_contigs"

# Display available GPUs
GPUS_AVAILABLE=$(nvidia-smi --query-gpu=index --format=csv,noheader | tr '\n' ' ')
echo "GPUs available: $GPUS_AVAILABLE"

# Calculate number of designs per GPU
GPU_n=$(echo $GPUS_AVAILABLE | wc -w)
designs_per_GPU=$((designs_n / GPU_n))

# Launch one job per GPU
for GPU_ID in $GPUS_AVAILABLE; do

    echo "Using $GPU_ID"
    (
        export CUDA_VISIBLE_DEVICES=$GPU_ID
        JOB_DIR="${SLURM_JOB_ID}_gpu${GPU_ID}"
        mkdir -p "$JOB_DIR"

        # Run partial diffusion
        echo "Running Partial Diffusion, generating in GPU ${GPU_ID}"
        /apps/rosetta/RFDifussion/scripts/run_inference.py inference.input_pdb="$input_pdb" inference.output_prefix="${JOB_DIR}/output/design" contigmap.contigs="$rfd_contigs" inference.num_designs="$designs_per_GPU" diffuser.partial_T="$noise_steps" denoiser.noise_scale_ca="$noise_scale" denoiser.noise_scale_frame="$noise_scale"  > "${JOB_DIR}/partial.out" 2> "${JOB_DIR}/partial.err"
        
        echo "Run ${JOB_DIR} finished"
    ) & 

done

wait
echo "All jobs completed"