#!/bin/bash

#Default option
pkg_manager="conda"


while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -p|--pkg) pkg_manager="$2" ; shift ;;
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
done

#Output selected package manager

echo -e "Package manager selected: $pkg_manager"

# set paths needed for installation and check for conda installation
install_dir=$(pwd)
CONDA_BASE=$(conda info --base 2>/dev/null) || { echo -e "Error: conda is not installed or cannot be initialised."; exit 1; }
echo -e "Conda is installed at: $CONDA_BASE"

### Start watcher install
echo -e "Installing watcher environment\n"
$pkg_manager create --name watcher python=3.10 -y || { echo -e "Error: Failed to create watcher conda environment"; exit 1; }
conda env list | grep -w 'watcher' >/dev/null 2>&1 || { echo -e "Error: Conda environment 'watcher' does not exist after creation."; exit 1; }

# Load newly created watcher environment
echo -e "Loading watcher environment\n"
source ${CONDA_BASE}/bin/activate ${CONDA_BASE}/envs/watcher || { echo -e "Error: Failed to activate the watcher environment."; exit 1; }
[ "$CONDA_DEFAULT_ENV" = "watcher" ] || { echo -e "Error: The watcher environment is not active."; exit 1; }
echo -e "watcher environment activated at ${CONDA_BASE}/envs/watcher"

# install required conda packages
echo -e "Instaling conda requirements\n"
$pkg_manager install pyrosetta --channel https://conda.graylab.jhu.edu -y  || { echo -e "Error: Failed to install conda packages."; exit 1; }
pip install CodonTransformer
pip install dash
pip install dash-bio
pip install biopython

# make sure all required packages were installed
required_packages=(pyrosetta codontransformer dash dash-bio biopython )
missing_packages=()

# Check each package
for pkg in "${required_packages[@]}"; do
    conda list "$pkg" | grep -w "$pkg" >/dev/null 2>&1 || missing_packages+=("$pkg")
done

# If any packages are missing, output error and exit
if [ ${#missing_packages[@]} -ne 0 ]; then
    echo -e "Error: The following packages are missing from the environment:"
    for pkg in "${missing_packages[@]}"; do
        echo -e " - $pkg"
    done
    exit 1
fi

echo -e "Successfully finished watcher enviroment installation!\n"
echo -e "Activate environment using command: \"$pkg_manager activate watcher\""
