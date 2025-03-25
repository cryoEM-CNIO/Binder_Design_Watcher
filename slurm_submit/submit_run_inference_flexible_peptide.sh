#!/bin/bash

# Number of desired cpus:
#SBATCH --nodes=1
#SBATCH -p cchacon
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
        --output_prefix)
            output_prefix="$2"
            echo "$2"
            shift # Shift past the argument value
            ;;
        --input_pdb)
            input_pdb="$2"
            echo "$2"
            shift
            ;;
        --contigmap_descriptor)
            contigmap_descriptor="$2"
            echo "$2"
            shift
            ;;
        --hotspots_descriptor)
            hotspots_descriptor="$2"
            echo "$2"
            shift
            ;;
        --designs_n)
            designs_n="$2"
            shift
            ;;
        --ckp)
            ckp="$2"
            shift
            ;;
        --conformation)
            conformation="$2"
            shift
            ;;
        --mask)
            mask="$2"
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
machine=`hostname`
echo "Current machine $machine"
echo "You chose the following input values:"
echo "  - output prefix: $output_prefix"
echo "  - input PDB: $input_pdb"
echo "  - contigmap descriptor: $contigmap_descriptor"
echo "  - hotspots descriptor: $hotspots_descriptor"
echo "  - number of designs: $designs_n"
echo "  - checkpoint used: $ckp"
echo "  - Mask conformation: $conformation"
echo "  - Structure Mask: $mask"

# Activate environment
source /apps/profile.d/load_all.sh
conda activate SE3nv

# Run!
if [ -z "${hotspots_descriptor+x}" ]; then
    /apps/rosetta/RFDifussion/scripts/run_inference.py inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor" inference.num_designs=$designs_n denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 inference.ckpt_override_path="$ckp"
else
    if [ ! -z "${mask}" ]; then
        if [ "$conformation" == "A" ]; then
            echo "Flexible peptide path selected: Alpha helix conformation "
            /apps/rosetta/RFDifussion/scripts/run_inference.py inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor" ppi.hotspot_res="$hotspots_descriptor" inference.num_designs=$designs_n denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 inference.ckpt_override_path="$ckp" +contigmap.inpaint_str="$mask" scaffoldguided.scaffoldguided=True +contigmap.inpaint_str_helix="$mask" 
        elif [ "$conformation" == "B" ]; then
        echo "Flexible peptide path selected: Beta Strand conformation "
            /apps/rosetta/RFDifussion/scripts/run_inference.py inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor" ppi.hotspot_res="$hotspots_descriptor" inference.num_designs=$designs_n denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 inference.ckpt_override_path="$ckp" +contigmap.inpaint_str="$mask" scaffoldguided.scaffoldguided=True +contigmap.inpaint_str_strand="$mask" 
        else
            echo "Flexible peptide path selected: No Given conformation "
            /apps/rosetta/RFDifussion/scripts/run_inference.py inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor" ppi.hotspot_res="$hotspots_descriptor" inference.num_designs=$designs_n denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 inference.ckpt_override_path="$ckp" +contigmap.inpaint_str="$mask"
        fi 
    else
        echo "Normal binder generation"
        /apps/rosetta/RFDifussion/scripts/run_inference.py inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor" ppi.hotspot_res="$hotspots_descriptor" inference.num_designs=$designs_n denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 inference.ckpt_override_path="$ckp"
    fi
fi

echo "done"
