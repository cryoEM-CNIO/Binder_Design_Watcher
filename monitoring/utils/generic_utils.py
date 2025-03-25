'''
This is where all the functions of the watcher are going to be stored in order to make the file more readable

'''
import pandas as pd
import os
import re
import time 
import numpy as np
import glob 
import argparse
import subprocess
import plotly.express as px
import plotly.graph_objects as go
from Bio import SeqIO

#Count the residues in chain in order to get the length
#Oye, muy bien pensado al que haya escrito esto, chapeau
def count_residues_in_chain(pdb_file, chain_id='A'):
    '''
    Function to compute the lenght of the binders. Since the backbone is composed of Gly, the number of residues is computed as the number of atoms divided by 4 (the number of atoms in a _Gly residue) (chapeau Nayim cause this is faster than the old function and very witty).

    Input:

    pdb_file --> Path to the pdb file whose length is being computed 

    chain_id --> Chain ID whose length we want to compute, always pointing at A which is the binder chain in RFD

    Output:

    residue_count/4 --> The number of residues of the binder  
    '''
    residue_count = 0
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line.split()[4] == chain_id:
                residue_count += 1
    return residue_count/4


def trim_substituted(text):
    '''
    Trim descriptions of silents to get original pdb designs names. This code is used to get the original design label. This is specially useful for the ngl visualization

    Input:

    text --> The description in the merged_df

    Output:

    whole_path --> Path to the original pdb in the output
    '''
    file = text.split('_substituted')[0]+'.pdb'
    whole_path = 'output/'+file.split('_design')[0]+'/'+file
    return whole_path

def merge_csv_files(directory, input_pdb_path):
    '''
    Function to merge all the scoring data (the one from AF2IG and from scoring.py into one df)

    Input:

    directory --> Directory where the data is stored (before the output)
    
    input_pdb_path --> Path to the input PDB is stored (For comparisons in scaffolding and partial diffusion cases, important because )
    
    Output:

    merged_df --> A df that integrates all the metrics from AF2IG and PyRosetta

    '''

    ################################################
    # SC reading no longer needed, scoring does it #
    # Kept for old projects, remove it in producti #
    ################################################
    try:
        input_pdb_path=f'{directory}/{input_pdb_path}'
        df_whole=pd.read_csv(f'{directory}/Scoring_Stats.csv')
        if 'time' not in df_whole.columns:
            df_list=[]
            working_directory=f'{directory}/output/'

            for root, dirs, files in os.walk(working_directory):
                for file in files:
                    if file.endswith('.sc'):
                        file_path = os.path.join(root, file)
                        df = pd.read_table(file_path, sep=r'\s+', encoding='utf-8')
                        df_list.append(df)
            if df_list:
                merged_df = pd.concat(df_list, ignore_index=True)
                merged_df=pd.merge(merged_df, df_whole, on='description')
                merged_df['original_design'] = merged_df['description'].apply(trim_substituted)
                input_df=pd.DataFrame(get_input_data(input_pdb_path))
                merged_df=pd.concat([merged_df,input_df], ignore_index=True)
                return merged_df
            else:
                return pd.DataFrame(columns=['plddt_binder','pae_interaction', 'CUTRE','dG','dSASA', 'Shape_complementarity', 'Packstat', 'dG_SASA_ratio', 'SAP', 'binder_int_hyd', 'binder_surf_hyd', 'interface_hbonds', 'interface_unsat_hbonds' ])
        else:
            df_whole['original_design'] = df_whole['description'].apply(trim_substituted)
            input_df=pd.DataFrame(get_input_data(input_pdb_path))
            merged_df=pd.concat([df_whole,input_df], ignore_index=True)
            return merged_df
    except FileNotFoundError:
         return pd.DataFrame(columns=['plddt_binder','pae_interaction', 'CUTRE','dG','dSASA', 'Shape_complementarity', 'Packstat', 'dG_SASA_ratio', 'SAP', 'binder_int_hyd', 'binder_surf_hyd', 'interface_hbonds', 'interface_unsat_hbonds' ])

# Function to track job status
def track_job_status(directory):
    '''
    This function tracks a job status so can be checked in a table below the scatter plot

    The soluble pmpnn and the scoring probably can be added 

    Input:
    
    directory --> It takes as input the working directory, including the output

    Output:

    status_df --> Returns a df with the status of each run, which can be RFDm pMPNN, AF2, waiting or failed  
    '''

    job_status = []
    for root, dirs, files in os.walk(directory):
        for dir_name in dirs:
            job_path = os.path.join(root, dir_name)
            trj_path = os.path.join(job_path, "trj")  # Path to the "trj" directory

            job_record = {"job": dir_name, "status": "Waiting"}  # Default status

            # Check for RFD status
            if os.path.exists(trj_path) and os.path.isdir(trj_path):
                job_record["status"] = "RFD" 

            # Check for pMPNN status
            input_files = [f for f in os.listdir(job_path) if f.endswith('_input.silent.idx')]
            if input_files:
                job_record["status"] = "pMPNN"

            # Check for AF2 status
            af2_files = [f for f in os.listdir(job_path) if f.endswith('_out.silent.idx')]
            if af2_files:
                job_record["status"] = "AF2"

            # Check if the job has finished
            done_file = [f for f in os.listdir(job_path) if f == dir_name + '_done']
            sc_files = [f for f in os.listdir(job_path) if f.endswith('.sc')]

            if done_file:
                if sc_files:  # If there are *.sc files, consider it finished
                    job_record["status"] = "Finished"
                else:  # If there are no *.sc files, mark as FAILED
                    job_record["status"] = "FAILED"

            job_status.append(job_record)
        break  # Exit after the first iteration to limit depth
    
    for record in job_status:
        job_name = record["job"]
        numeric_part = int(re.search(r'\d+', job_name).group())  # Extract numeric part
        record["sort_key"] = numeric_part  # Add as sort_key

    # Sort by this new key before passing to DataTable
    status_df = pd.DataFrame(job_status).sort_values(by="sort_key")
    return status_df

#Function to get the data from the input
def get_input_data(input_pdb_path):
    '''
    Reads the last lines of the pdb used as input and stores the information in a df. The input must be the output from the watcher (to have all the metrics)
    or the output of a previous RFD run (only has the "soft filter" metrics)
    This is meant to be used in the case of scaffold or partial diffusion

    Input:

    input_pdb_path --> The path of where the input is stored

    Output:

    binding_analysis_dict --> A dictionary containing all the metrics of the initial input
    '''

    #Get the file working with regular expression using glob
    try:
        input_pdb_path=glob.glob(input_pdb_path)[0]
    except IndexError:
        input_pdb_path=input_pdb_path
    # Load dictionary
    binding_analysis_dict = {
        'pae_interaction':[],
        'pae_binder':[],
        'pae_target':[],
        'plddt_total':[],
        'plddt_binder':[],
        'plddt_target':[],
        'description': [],
        'original_design': [],
        'CUTRE': [],
        'dG': [],
        'dSASA': [],
        'Shape_complementarity': [],
        'Packstat': [],
        'dG_SASA_ratio': [],
        'SAP': [],
        'binder_int_hyd': [],
        'binder_surf_hyd': [],
        'interface_hbonds': [],
        'interface_unsat_hbonds': []
    }
    try:
        with open(f'{input_pdb_path}', 'r') as file:
            lines=file.readlines()
            for key in binding_analysis_dict.keys():
                for line in lines:
                # Iterate over the dictionary keys
                    if line.startswith(key):
                        # Extract the value after the key and any whitespace
                        value = line[len(key):].strip()
                        # Append the value to the corresponding key's list
                        binding_analysis_dict[key].append(float(value))
                        break
        
        binding_analysis_dict['original_design'].append(input_pdb_path)
        binding_analysis_dict['description'].append(input_pdb_path.split('/')[1].split('.')[0])
        #check if some of the keys of the dictionary are empty and add a nan
        for key, value in binding_analysis_dict.items():
            if isinstance(value, list) and not value:  # Check if it's an empty list
                value.append(np.nan)
        print("Data loaded successfully.")
    except FileNotFoundError:
        print(f"File {input_pdb_path} not found.")
        for key in binding_analysis_dict.keys():
            binding_analysis_dict[key].append(np.nan)

    return binding_analysis_dict


#Filtering of DataFrame to get only hits according to the metrics
def filtering_df(merged_df,pae_interaction_thres, CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres):
    '''
    --> Function to filter the df for the selection between the threshold. Separated to make code more legible

    Input:
    
    merged_df --> The dataframe with all the data
    All the thresholds --> (list of two elements containing the maximum and the minimum thresholds) 

    Output:

    filtered_df --> A df containing only the designs that fulfill all the metrics 
    
    '''

    filtered_df=merged_df[merged_df['pae_interaction'].between(pae_interaction_thres[0], pae_interaction_thres[1])]
    filtered_df=filtered_df[filtered_df['CUTRE'].between(CUTRE_thres[0], CUTRE_thres[1])]
    filtered_df=filtered_df[filtered_df['plddt_binder'].between(plddt_binder_thres[0], plddt_binder_thres[1])]
    filtered_df=filtered_df[filtered_df['Shape_complementarity'].between(shape_complementarity_thres[0], shape_complementarity_thres[1])]
    filtered_df=filtered_df[filtered_df['dSASA'].between(dsasa_thres[0], dsasa_thres[1])]
    filtered_df=filtered_df[filtered_df['interface_hbonds'].between(interface_hbond_thres[0], interface_hbond_thres[1])]
    filtered_df=filtered_df[filtered_df['interface_unsat_hbonds'].between(interface_unsat_hbond_thres[0], interface_unsat_hbond_thres[1])]
    filtered_df=filtered_df[filtered_df['binder_surf_hyd'].between(binder_surf_hyd_thres[0], binder_surf_hyd_thres[1])]

    return filtered_df



def get_working_directories(working_dir):
    """
    Get a list of directories within the parent directory that contain a folder named 'output'.
    
    Input:
    parent_dir (str): Path to the parent directory.
    
    Output:
    list: A list of paths to directories containing a folder named 'output'.
    """
    directories_list = []

    for root, dirs, files in os.walk(working_dir):
        if 'Scoring_Stats.csv' in files:  # Check if 'Scoring_Stats.csv' is one of the directories, since this file is needed for the plotting
            directories_list.append(root)

    return directories_list












