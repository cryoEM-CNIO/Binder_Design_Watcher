#Manually set the paths for every other software you have installed and environment name

RFD_PATH="/path/to/RFDiffusion"
PMPNN_PATH="/path/to/dl_binder_design" #If you are usign the tools from nrbennett repo, it should be the same as SILENT and AF2IG
AF2IG_PATH="/path/to/dl_binder_design" #If you are usign the tools from nrbennett repo, it should be the same as PMPNN and SILENT
MICRORUN_PATH="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
SILENT_PATH="/path/to/dl_binder_design" #If you are usign the tools from nrbennett repo, it should be the same as PMPNN and AF2IG

# MANUALLY SET ENVIRONMENTS NAMES (example names set)

RFD_ENV="SE3nv"
PMPNN_ENV="dl_binder_design"
AF2_ENV="dl_binder_design"
MICRORUN_ENV="watcher"


#MANUALLY SET SBATCH CONFIGURATIONS (examples in place)

NODES=1
PARTITION=Binder_design
CPUS_PER_GPU=12
GRES=gpu:1