#!/usr/bin/env python3
'''
20241106

Code to compute all by all rmsd and clusterize them into X different clusters using hierarchichal clustering

Input:

--directory: Directory where the pdbs are (eg the hits folder)
-- t: Number of clusters

Output:

-- f'{args.directory}/means_by_cluster.csv' : A df with the metrics of each cluster
-- f'{args.directory}/hits_clustered.csv' : The df we all know with a column called cluster for the cluters
'''




from pyrosetta import *
import glob 
import re
import os
import pandas as pd
import numpy as np
import pyrosetta.toolbox
from superimpose import superpose_pose_by_chain #This is a function for another code, recommended to have this code in the same folder as the one you are reading
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('--directory', '-d', help='Directory where the pdb files are stored')
parser.add_argument('--t', '-t', help='number of clusters to make', type=int, default=4)
args=parser.parse_args()

#Initiate pyrosetta
pyrosetta.init()

#Get the files
file_list=glob.glob(f'{args.directory}/*pdb')
print(file_list)

#Get the numbering 
number_of_files=len(file_list)
print(f'TOTAL NUMBER OF FILES: {number_of_files}\n')
print(f'TOTAL NUMBER OF CALCULATIONS: {number_of_files*(number_of_files-1)/2}\n')
#Create the nxn matrix
rmsd_matrix=pd.DataFrame(np.zeros((number_of_files, number_of_files)), index=file_list, columns=file_list)

#Pattern

pattern=r'.*/(run_\d*_design_\d*.*af2pred).*'


#Compute and assign the pairwise RMSD to the matrix. BEWARE! There is not very well explained doc about rosetta online, so not 100% sure the RMSD it gives you is the tru RMSD
for i in range(number_of_files):
    structure1=file_list[i]
    pose1=pose_from_pdb(structure1) 
    for j in range(i+1,number_of_files):
        structure2=file_list[j]
        pose2=pose_from_pdb(structure2)
        superpose_pose_by_chain(pose1, pose2, chain='B') #Aligning
        pose1_binder=pose1.split_by_chain(1)
        pose2_binder=pose2.split_by_chain(1)
        rmsd_CA_fs=pyrosetta.rosetta.core.scoring.CA_rmsd(pose1,pose2)
        rmsd_CA_A=pyrosetta.rosetta.core.scoring.CA_rmsd(pose1_binder,pose2_binder) #Compute
        rmsd_matrix.iloc[i,j]=rmsd_CA_A
        rmsd_matrix.iloc[j,i]=rmsd_CA_A

np.fill_diagonal(rmsd_matrix.values, 0)

rmsd_matrix.to_csv(f'{args.directory}/RMSD_matrix.csv')
# Convert RMSD matrix to condensed form (1D array)
rmsd_condensed = squareform(rmsd_matrix)

# Perform hierarchical clustering
Z = linkage(rmsd_condensed, method='average')  # Use 'average' or 'complete' linkage

# Assign clusters (4 clusters in this case)
cluster_labels = fcluster(Z, t=args.t, criterion='maxclust')

#Create a dictionary to store the names and the cluster they belong

dictionary={
    'description':[],
    'cluster':[]
}


for i in range(number_of_files):
    print(re.search(pattern, file_list[i]).group(1))
    dictionary['description'].append(re.search(pattern, file_list[i]).group(1))
    dictionary['cluster'].append(cluster_labels[i])

print(dictionary)
cluster_df=pd.DataFrame(dictionary)
hits_df=pd.read_csv('hits.csv')

merged_df=pd.merge(cluster_df,hits_df, on='description')

# merged_df=merged_df.drop(['SCORE'], axis=1)
# merged_df=merged_df.drop(['Unnamed: 0'], axis=1)


CUTRE_mean_by_cluster=merged_df.groupby('cluster')['CUTRE'].mean()
PAE_mean_by_cluster=merged_df.groupby('cluster')['pae_interaction'].mean()
PLDDT_mean_by_cluster=merged_df.groupby('cluster')['plddt_binder'].mean()
length_mean_by_cluster=merged_df.groupby('cluster')['length'].mean()
proteins_by_cluster=merged_df.groupby('cluster').size().rename('Number_of_Proteins')

means_df=CUTRE_mean_by_cluster.to_frame().join([PAE_mean_by_cluster, PLDDT_mean_by_cluster, length_mean_by_cluster, proteins_by_cluster])
means_df.to_csv(f'{args.directory}/means_by_cluster.csv')
merged_df.to_csv(f'{args.directory}/hits_clustered.csv')

for cluster in merged_df['cluster']:
    with open(f'cluster{cluster}.txt', 'w') as file:
        for protein in merged_df['description'][merged_df['cluster']==cluster]:
            file.write(f'{protein}.pdb\n')


print('ALL RMSD COMPUTED: DONE!\n')
print('\n')
print('METRICS PER CLUSTER SAVED AT means_by_cluster.csv\n')
print('HITS AND CLUSTER SAVED AT hits_clustered.csv\n')
print('PROTEIN PDBS PER CLUSTER SAVED AT clusterX.txt')
