#!/bin/bash


# Parse command-line arguments  
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --number) number="$2" ; shift ;;
        --run) run="$2" ; shift ;;
        --directory) directory="$2" ; shift ;; 
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift # Shift past the current argument
done

#Load all variables

source $directory/config.sh
conda activate $MICRORUN_ENV



# Display the parsed values
machine=`hostname`
echo "Current machine $machine"

output="output/run_${run}/run_${run}_done"

echo "input_silent: $output"

# Run
touch $output
python3 $MICRORUN_PATH/microrun/scripts/run_ending.py --number "$number"

echo "Touched $output"
