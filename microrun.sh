#! /bin/bash

#Get script dir and load all the variables
SCRIPT_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
source $SCRIPT_DIR/config.sh 

#Checking everything
echo $SCRIPT_DIR
echo $MICRORUN_PATH
echo $RFD_PATH
echo $PMPNN_PATH


# Master script for deep searches in microruns

# 1: Get all info needed
##Set defaults
partial_diff="False"
noise_steps=20
pmp_nseqs=1
rfd_ndesigns=10
pmp_relax_cycles=1
noise_scale=1
ckp="$RFD_PATH/models/Complex_base_ckpt.pt"
node=''
soluble="False"
distance=10
hits_number=100
core=0.05
residues=None


## Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input) input="$2" ; shift ;;
        -t|--template) template="$2" ; shift ;;
        -m|--max_threads) max="$2" ; shift  ;;    
        -c|--rfd_contigs) rfd_contigs="$2" ; shift  ;;    
        -h|--rfd_hotspots) rfd_hotspots="$2" ; shift  ;;    
        -np|--pmp_nseqs) pmp_nseqs="$2" ; shift  ;;    
        -rc|--pmp_relax_cycles) pmp_relax_cycles="$2" ; shift  ;;   
        -pd|--partial_diff) partial_diff="$2" ; shift  ;; 
        -nst|--noise_steps) noise_steps="$2" ; shift  ;;
        -nsc|--noise_scale) noise_scale="$2" ; shift  ;;
        -ck|--ckp) ckp="$2" ; shift  ;; #Add the path to the checkpoint to add weight toward some fold
        -w|--node) node="$2" ; shift  ;; # Provide a specific name of a node to submit to this node with -w. If not provided it will be as usual.
        -sp|--soluble_pMPNN) soluble="$2" ; shift  ;; #Add a pMPNNsol step to the non interacting interface to increase the solubility
        -d|--distance) distance="$2" ; shift  ;; # Distance to fix residues in the interaface
        -hn|--hits_number) hits_number="$2" ; shift ;;
        -cr|--core) core="$2" ; shift ;; # Proportion of core residues to make the filtering
        -re|--residues) residues="$2" ; shift ;; #Residues index to fix, useful for scaffolding
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift  # Shift past the current argument
done


mkdir -p ./slurm_logs
mkdir -p ./output

last_run_folder=$(ls -d "./output/run_"* 2>/dev/null | sort -V | tail -n 1 | sed 's#./output/run_##')

if [[ -n "$last_run_folder" ]]; then
    i="$last_run_folder"
else
    i=0
fi

## Contigs Getter, automatic for partial diffusion

if [ $partial_diff = "True" ]; then 
        rfd_contigs=$(python3 $MICRORUN_PATH/scripts/contigs_map_getter.py --input "$input" --partial_diff "$partial_diff")
        echo "$rfd_contigs"
fi

## Prepare Folder & Variables
echo "Preparing JSON to save the run variables"
python3 $MICRORUN_PATH/microrun/scripts/json_variable_generation.py --input "$input" --template "$template" --max_threads "$max" --rfd_contigs "$rfd_contigs" --rfd_hotspots "$rfd_hotspots" --rfd_ndesigns "$rfd_ndesigns" --pmp_nseqs "$pmp_nseqs" --pmp_relax_cycles "$pmp_relax_cycles" --partial_diff "$partial_diff" --noise_steps "$noise_steps" --noise_scale "$noise_scale" --ckp "$ckp" --soluble_pMPNN "$soluble" --distance "$distance" --core "$core" --residues "$residues"

old_i=1

# RUN

while [ ! -f 'campaign_done' ]; do
    i=$((i+1))
    echo "Starting cycle $i"
    mkdir -p ./output/run_$i

    # Set a default value for previous if i is smaller than max (initial cycles)
    if [ "$i" -le $max ]; then
        previous=0
    else
        previous=$((i - $max))
    fi
    waitfor="output/run_${previous}/run_${previous}_done"
    # Skip waiting when previous is 0
    if [ "$previous" -ne 0 ]; then
        # Wait until the condition is met
        while [ ! -e "$waitfor" ]; do
            echo "Waiting for previous jobs to complete: $waitfor"
            sleep 60
        done
    # Now that the condition is met or previous is 0, proceed to the following code
    fi  

sbatch -w "$node" --nodes="$NODES" -p "$PARTITION" --open-mode=append --gres="$GRES" --exclusive --cpus-per-gpu="$CPUS_PER_GPU" -o slurm_logs/%j.out -e slurm_logs/%j.err \
       "$MICRORUN_PATH/microrun/slurm_submit/submit_master.sh" --input "$input" --template "$template" --run "$i" --rfd_contigs "$rfd_contigs" --rfd_ndesigns "$rfd_ndesigns" \
       --pmp_nseqs "$pmp_nseqs" --pmp_relax_cycles "$pmp_relax_cycles" --partial_diff "$partial_diff" --noise_steps "$noise_steps" --noise_scale "$noise_scale" --ckp "$ckp" \
       --core "$core" --residues "$residues" --hits_number "$hits_number" --rfd_hotspots "$rfd_hotspots" --directory "$SCRIPT_DIR"

done

echo "Campaign finshed"
