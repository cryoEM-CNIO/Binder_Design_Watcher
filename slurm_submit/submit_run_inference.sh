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

export HYDRA_FULL_ERROR=1

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --output_prefix)output_prefix="$2" ; shift ;;
        --input_pdb) input_pdb="$2" ; shift ;;
        --contigmap_descriptor) contigmap_descriptor="$2" ; shift ;;
        --hotspots_descriptor) hotspots_descriptor="$2" ; shift ;;
        --designs_n) designs_n="$2" ; shift ;;
        --ckp) ckp="$2" ; shift ;;
        --partial_diff) partial_diff="$2" ; shift ;;
        --noise_steps) noise_steps="$2" ; shift ;;
        --noise_scale) noise_scale="$2" ; shift ;;
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift 
done

# Display the parsed values
machine=`hostname`
echo "Current machine $machine"
echo "You chose the following input values:"
echo "  - output prefix: $output_prefix"
echo "  - input PDB: $input_pdb"
echo "  - contigmap descriptor: $contigmap_descriptor"
echo "  - hotspots descriptor: $hotspots_descriptor"
echo "  - number of designs: $designs_n"
echo "  - checkpoint used: $ckp"
echo "  - Partial Diffusion: $partial_diff"
echo "  - noise scale (only if PD): $noise_scale"
echo "  - noise_steps (only if PD): $noise_steps"

# Activate environment
source /apps/profile.d/load_all.sh
conda activate SE3nv

# Run!
if [ $partial_diff = "True" ]; then
    /apps/rosetta/RFDifussion/scripts/run_inference.py inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor"  inference.num_designs="$designs_n" diffuser.partial_T="$noise_steps" denoiser.noise_scale_ca="$noise_scale" denoiser.noise_scale_frame="$noise_scale"
else
    if [ -z "${hotspots_descriptor+x}" ]; then
        /apps/rosetta/RFDifussion/scripts/run_inference.py inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor" inference.num_designs=$designs_n denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 inference.ckpt_override_path="$ckp"
    else
        /apps/rosetta/RFDifussion/scripts/run_inference.py inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor" ppi.hotspot_res="$hotspots_descriptor" inference.num_designs=$designs_n denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 inference.ckpt_override_path="$ckp"
    fi
fi
echo "done"
