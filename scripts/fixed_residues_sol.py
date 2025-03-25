#!/usr/bin/env python3

from pyrosetta import *
import argparse
import pandas as pd
import re
import json
import numpy as np
import glob 


init()

# Protein name
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='input file path')
parser.add_argument('--distance', type=float, default=10, help='distance between atoms to consider a residue is interacting with the target. The default is 10 ')
args = parser.parse_args()
#Load files
file=args.input
distance=args.distance

#Load list
close_residues_binder_dict={}

# CLOSE RESIDUES SELECTION 

#Variable and classes loading

poses=poses_from_silent(file)

for pose in poses:
    NRS_binder=pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    NRS_target=pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    NRS_target.set_distance(10)
    NRS_binder.set_distance(10)
    interface='A_B'

    #Save the protein name in the dict
    protein_name=pose.pdb_info().name()
    close_residues_binder_dict[protein_name]=[]

    #Get the binder start and ending numeration
    start_residue_binder=pose.conformation().chain_begin(1)
    end_residue_binder=pose.conformation().chain_end(1)

    #Set focus and search neighbors of the target
    NRS_binder.set_focus(f'{start_residue_binder}-{end_residue_binder}')
    neighbors_binder=NRS_binder.apply(pose)



    #Get the target numeration
    start_residue_target=pose.conformation().chain_begin(2)
    end_residue_target=pose.conformation().chain_end(2)

    #Set focus and search
    NRS_target.set_focus(f'{start_residue_target}-{end_residue_target}')
    neighbors_target=NRS_target.apply(pose)

    #Save the neighbors of the binder
    close_residues_binder=[]
    for j in range(start_residue_binder, end_residue_binder):
        if neighbors_target[j]:
            close_residues_binder_dict[protein_name].append(j)
    

    #########################################################################

print(f"Creating JSON file of fixed residues")
directory = os.path.dirname(file)

file_path = os.path.join(directory, "fixed_residues.json")

# Writing the list to a JSON file
with open(file_path, 'w') as json_file:
    json.dump(close_residues_binder_dict, json_file, indent=4)  # 'indent' makes the file more readable
print(f"List of numbers saved to {file_path}")

