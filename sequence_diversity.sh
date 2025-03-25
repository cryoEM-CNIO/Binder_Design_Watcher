#!/bin/bash

# Load defaults 
fr=1
nseqs=1
fixed="None"

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --input) input="$2" ; shift ;;
        --threads) threads="$2" ; shift ;;
        --max) max="$2" ; shift ;;
        --fixed) fixed="$2" ; shift ;;
        --nseqs) nseqs="$2" ; shift ;;
        --fr) fr="$2" ; shift ;;
        *)
        esac
        shift
done

mkdir -p "output"
mkdir -p "slurm_logs"

last_run_folder=$(ls -d "./output/run_"* 2>/dev/null | sort -V | tail -n 1 | sed 's#./output/run_##')

if [[ -n "$last_run_folder" ]]; then
    i="$last_run_folder"
else
    i=0
fi

## Fix residues

if [ $fixed != "None" ]; then 
    python3 /apps/scripts/protein_design/scripts/fixing_residues.py --fixed "$fixed" --pdb_input "$input"
fi

## create silent file
/apps/rosetta/dl_binder_design/include/silent_tools/silentfrompdbs  "$input" > "initial_input.silent"


#while loop to generate diversity
while true;do
    i=$((i+1))
    if [ "$i" -lt "$threads" ]; then
        previous=0
    else
        previous=$((i - $threads))
    fi
    
    # Variable generation and definition

    mkdir -p "output/run_${i}"

    waitfor="output/run_${previous}/run_${previous}_done"

    if [ "$previous" -ne 0 ]; then 
        
        while [ ! -e "$waitfor" ]; do 
            echo "Waiting for previous job to complete: $waitfor"
            sleep 60
        done
    
    fi

    sbatch /apps/scripts/protein_design/slurm_submit/submit_sequence_diversity.sh --run "$i" --nseqs "$nseqs" --fr "$fr" --fixed "$fixed"

    total_seqs_generated=$(($i*4*$nseqs))
    if [ $total_seqs_generated -gt $max ]; then 
        break 
    fi  
done
counter=$((counter+1))

echo " #######################"
echo "- DIVERSITY GENERATION CONCLUDED"
echo "- ${total_seqs_generated} sequences generated in total"
echo " #######################"