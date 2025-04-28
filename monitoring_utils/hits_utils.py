import pandas as pd
import os
import re
import warnings
import time 
import numpy as np
import glob 
import argparse
import json
import subprocess
import random
import math
import plotly.express as px
import plotly.graph_objects as go
import warnings
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# SeqIO raises some warnings that not affect the outcome, this silence them
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# Function to get best hit name
def get_hit_names(filtered_df, xaxis_value):
    '''
    Function to get all the hit names. This is employ to get the hit names for the ngl representation dropdown
    Probably could be fused with the enxt function

    Input: 

    filtered_df --> DF filtered with only those designs that fulfill all the filters

    xaxis_value --> Variable which is employ to order the designs in the dropdown
    
    Output:

    hit_names --> A dictionary of the hit names, ordered in ascending order using the x variable
    '''
    hit_names=[]
    if not filtered_df.empty:
        sorted_df = filtered_df.sort_values(by=xaxis_value, ascending=True)
        for index, row in sorted_df.iterrows():
            hit_names.append(row['description'])
        return hit_names
    else:
        return ["No hits found under the specified conditions"]
    
# Function to get best hit path
def get_design_file_path_and_name(hits_names, directory,input_pdb_path):
    '''
    Gets a design file path and name for its representaion in ngl

    Input:

    description --> identifier of the hit name

    directory --> directory where we are working /path/to/output/.
    
    Output:

    Returns the data path and its identification for the ngl representation  

    '''
    
    description = hits_names
    match = re.match(r"(run_\d+)_design_(\d+_substituted).*", description)
    match2 = re.match(r"Input", description)
    if not match and not match2:
        print("Invalid description format")
        # Returning None for both values if the format is invalid
        return None, None
    elif match2:
        input_path=os.path.join(directory, input_pdb_path)
        # If the description is "Input", return the directory and a placeholder filename
        data_path ='/'.join(glob.glob(input_path)[0].split('/')[:-1])
        filename = glob.glob(input_path)[0].split('/')[-1]
        return data_path, filename
    run_part, design_part = match.groups()
    # Constructing the file path using 'directory' which is now an absolute path
    working_directory = directory+'/output'
    data_path = os.path.join(working_directory, run_part)
    filename = f"{run_part}_design_{design_part}"

    return data_path, filename

#Function to add the stats of the hits to the PDB, as they are added in after AF2IG
#Makes PD comparison more comfortable
def add_stats_to_pdb(description,directory):
    '''
    add the stats to the pdb, in the same fashion af2IG does (adding the name and the metric separated by a whitespace)

    Input:

    description --> identifier of the hit to record the metrics inside them

    directory --> directory where we are working and were the Scoring_Stats.csv is

    Output:

    Writes the data at the end of the pdb file    
    '''

    df_rosetta = pd.read_csv(f'{directory}/Scoring_Stats.csv')
    design_metrics = df_rosetta[df_rosetta['description'] == description]
    design_metrics = design_metrics.drop(['close_residues_target', 'close_residues_binder'], axis=1)
    #The fastas has an abbreviated form of the name, so you have to crop it
    description_short=re.search(r'(run_\d+_design_\d+.*_dldesign_\d+).*', description).group(1)
    pdb_path = f'{directory}/hits/{description}.pdb'
    fasta_path=f'{directory}/hits/fastas/{description_short}.fasta'
    mw, ip, extinction_coefficient=param_stats(fasta_path)
    with open(pdb_path, 'a') as file:
        for column in design_metrics.columns:
            if column != 'description':  # Skip writing the 'description' column itself
                value = design_metrics[column].iloc[0]  # Extract the first value from the series
                file.write(f'{column} {round(value,4)}\n')
        file.write(f'molecular_weight {round(mw,4)}\n')
        file.write(f'isoelectric_point {round(ip,4)}\n')
        file.write(f'extinction_coefficient {round(extinction_coefficient,4)}\n')

#Function to extract the fasta seqs of the hits 
def extract_fasta_seq(description):
    '''
    Extract fasta sequences for order

    Input:

    description --> Hit name

    Output:

    {description}.fasta --> Fasta sequence of the hit
    '''

    # Define variables
    input_file=f'{description}.pdb'
    pattern = r'.*(run_\d+_design_\d+.*_dldesign_\d+)'
    directory = os.path.dirname(input_file)

    #Get fasta names 
    try:
        fasta_name = re.search(pattern, input_file).group(1)
    except AttributeError:
        fasta_name = input_file[:-4]
    
    
    #Extract sequence 
    with open(input_file, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            chain = record.annotations['chain']
            if chain == 'A':
                sequence = record.seq

    #write things down

    fasta_file = f'{directory}/fastas/{fasta_name}.fasta'
    with open(fasta_file, 'w') as fasta:
        fasta.write(f'>{fasta_name}\n{sequence}\n')

    return fasta_file
# Function to extract the DNA seq with the desired characteristics
def extract_dna_seq(input, output,organism, met, overhang_5 , overhang_3, length, GC,enzyme  ):
    '''
    Extract optimized dna seqs for order

    Input:

    input --> Fasta file with the aminoacidic sequence

    output --> Folder in which the dna seqs are gonna be stored

    organism --> Organism desired for codon optimization

    met --> Whether to add a methionine at the start of the sequence or not (default: True)

    overhang_5 --> Overhang sequence to add at the 5\' (default=Nothing) 

    overhang_3 --> Overhang sequence to add at the 3\' (default=Nothing)

    length --> Number of random bases that must be added to reach a certain size. Half of the bases are added at each extreme, between the overhangs and the sequence (default=0)        
    
    GC --> GC content of the random sequence, express in % (default=50)

    enzyme --> Restriction enzymes to check in the sequence
    Output:

    Nothing, it writes the DNA seq at the folder specified in the output key
    '''

    #Define a dictionary which is later gonna be used to store all the order information
    order_dictionary={
        'design_name':[],
        'organism':[],
        'met_added':[],
        '5\' overhang':[],
        '3\' overhang':[],
        'fasta':[],
        'dna_seq':[],
        'binder_order':[]
    }

    #Get the fasta sequence and design name
    with open (input, 'r') as fasta_file:
        for line in fasta_file.readlines():
            if not line.startswith('>'):
                sequence=line
            else:
                name=line 
                order_dictionary['design_name'].append(name[1:])
    
    #Add initial met 
    if met == True:
        sequence='M'+sequence
        met_fasta=input
        with open(met_fasta, 'w') as new_fasta_file:
            new_fasta_file.write(f'{name}')
            new_fasta_file.write(f'{sequence}')
        codon_transformer_input=met_fasta
    else:
        codon_transformer_input=input
            ##Variable preparation
    codon_transformer_output=os.path.join(output, codon_transformer_input.split('/')[-1].split('.')[0]+'_RevTrans.fasta')

    all_ok=False
    while all_ok == False:

        #run CodonTransformer

        command = (
            f"python3 /apps/scripts/protein_design/scripts/CodonTransformer_seq.py --protein '{sequence}' --organism '{organism}'"
        )

        # Execute the command
        process = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True
        )

            # Check for errors
        if process.returncode != 0:
            print("Error:", process.stderr)
            raise RuntimeError("CodonTransformer_process failed")

        dna_seq=process.stdout.strip().lower()


        ## Design random sequences for the 3 and 5 ends to reach the desired size
        length_to_add=int(length)-(len(dna_seq)+len(overhang_3)+len(overhang_5))
        if length_to_add > 0:
            print(f'A random sequence {length_to_add} residues long will be added')
        random_sequence_5,random_sequence_3= RandomSequenceGenerator(length_to_add, GC)


        full_seq=f'{random_sequence_5}'+f'{overhang_5}' + f'{dna_seq}'+ f'{overhang_3}' + f'{random_sequence_3}'
        full_seq,all_ok=check_enzyme_cut(enzyme,dna_seq)

    ## Open the RevTrans fasta_file

    with open (codon_transformer_output, 'w') as dna_seq_fasta:
        dna_seq_fasta.write(f'{name}')
        dna_seq_fasta.write(full_seq)
    
    
    print(f'Sequence RevTranslated and saved to: {codon_transformer_output}')
    print(f'The follwoing overhangs have been added\n' + f'5\' overhang: {overhang_5}\n' + f'3\' overhang: {overhang_3}')
    


def RandomSequenceGenerator(length, avoid_seqs,GC=50):
    '''
    Generate random DNA sequences with the desired GC proportion. The sequence is added symmetrically 
    at the 3' and 5' ends, with half the required length at each terminus.

    Input:
        length (int): Length of the random sequence to add. It is split equally between 5' and 3' ends.
        GC (float): Desired GC content in percentage (default=50).
        Avoid_seqs (list): List of sequences to avoid in the random sequence generation.

    Output:
        tuple: (sequence_5, sequence_3) Random sequences for the 5' and 3' ends.
    '''
    if length <= 0:
        return '', ''

    # Convert GC content percentage to probabilities
    GC_prob = GC / 100 / 2
    AT_prob = (1 - GC_prob * 2) / 2
    nucleotide_weights = [AT_prob, GC_prob, AT_prob, GC_prob]  # For 'a', 'g', 't', 'c'

    # Calculate the lengths for 5' and 3' ends
    length_5 = (length + 1) // 2  # Round up for the 5' end
    length_3 = length // 2       # Remaining for the 3' end
    
    # Initialize sequences
    seqs_avoided=False
    while seqs_avoided==False:
        # Initialize sequences
        sequence_5 = []
        sequence_3 = []

        print(f'Creating a random sequence {length} residues long')
        # Generate sequence for the 3' end
        for _ in range(length_3):
            sequence_3.append(random.choices(['a', 'g', 't', 'c'], weights=nucleotide_weights)[0])

        # Generate sequence for the 5' end
        for _ in range(length_5):
            sequence_5.append(random.choices(['a', 'g', 't', 'c'], weights=nucleotide_weights)[0])
        
        for seq in avoid_seqs:
            if seq.lower() in ''.join(sequence_5) or seq.lower() in ''.join(sequence_3):
                print(f'Sequence to avoid {seq} in the random sequence, starting again: :(')
                seqs_avoided=False
            else:
                seqs_avoided=True
                    
    # Convert list to string and return
    return ''.join(sequence_5), ''.join(sequence_3)


def generate_order_csv(extraction_list, hits_folder):
    '''
    Generates a csv which all the information required for ordering (in our case)

    Input:

    extraction_list --> List of designs whose structure and sequence is going to be extracted

    hits_folder --> Folder where all the hits information is being saved
    Output:

    order.csv -->  A csv file which details the path to the original structure, as well as to the aminoacidic and nucleotidic sequences, 
                   and the design identifier, aminoacidic sequence and nucleotidic sequence
    
    '''
    
        # Load dna_seqs
    data_dictionary={
        'design_path':[],
        'description':[],
        'fasta_path':[],
        'protein_sequence':[],
        'dna_path':[],
        'dna_sequence':[],
        'molecular_weight':[],
        'isoelectric_point':[],
        'molar_extinction_coefficient':[]
    }

    pattern = r'.*(run_\d+_design_\d+.*_dldesign_\d+).*'
    pattern_run=r'.*(run_\d+)_.*'


    #Load protein_seq
    for description in extraction_list:
        protein_path=f'{hits_folder}/fastas/{re.search(pattern,description).group(1)}.fasta' 
        with open(protein_path, 'r') as protein_file:
            for line in protein_file.readlines():
                if not line.startswith('>'):
                    protein_seq=line
        mw,ip,extinction_coefficient=param_stats(protein_path)
        #Load DNA seq
        dna_path=f'{hits_folder}/dna_seqs/{re.search(pattern,description).group(1)}_RevTrans.fasta'
        with open(dna_path, 'r') as dna_file:
            for line in dna_file.readlines():
                if not line.startswith('>'):
                    dna_seq=line

        run_number=re.search(pattern_run, description).group(1)
        path=f'output/{run_number}/{description}'
        data_dictionary['design_path'].append(path)
        data_dictionary['description'].append(re.search(pattern, description).group(1))
        data_dictionary['fasta_path'].append(protein_path)
        data_dictionary['protein_sequence'].append(protein_seq)
        data_dictionary['dna_path'].append(dna_path)
        data_dictionary['dna_sequence'].append(dna_seq)
        data_dictionary['molecular_weight'].append(mw)
        data_dictionary['isoelectric_point'].append(ip)
        data_dictionary['molar_extinction_coefficient'].append(extinction_coefficient)
    
    data_df=pd.DataFrame(data_dictionary)
    data_df.to_csv(f'{hits_folder}/order.csv')

def create_log_extraction(output_folder, extraction_list, met, organism, overhang_3, overhang_5, length, GC):
    '''
    Function to create a log file of the extraction. This file follows a JSON organization

    Input:

    output_folder --> Folder to store the log

    extraction_list --> List of designs to extract 
    
    met --> Whether if a Met has been added or not

    organism --> Organism for which the sequence has been optimized

    overhang_3 --> Sequence added to the 3\' overhang

    overhang_5 --> Sequence added to the 5\' overhang

    length --> Minimum length that the DNA seqs must reach (random seq added if the sequence length is inferior)

    GC -->   GC content of the random sequences added, in %

    Output:

    output_folder/extraction.log -->  A JSON file which details the conditions used for the extraction  

    '''
    extraction_log={
        'description_list':extraction_list, 
        'Methionine':[met],
        'Organism':[organism],
        '3_overhang':[overhang_3],
        '5_overhang':[overhang_5],
        'length':[length],
        'GC':[GC]
    }
    with open(f'{output_folder}/extraction.log', 'w') as outfile:
        json.dump(extraction_log,outfile)

def param_stats(fasta):
    '''
    This function is meant to check if there are aromatic residues in the protein to facilitate protein detection and quantification
    
    Input:

    fasta --> Fasta file with the protein sequence

    Output:

    mw --> molecular weight of the protein
    ip --> isoelectric point of the protein
    extinction_coefficient --> Extinction coefficient of the protein 

    '''

    #Extract the sequence

    with open(fasta, 'r') as file:
        for line in file.readlines():
            if not line.startswith('>'):
                sequence=line

    #open it in biopython suite
    X= ProteinAnalysis(sequence)

    #param calculation

    mw=X.molecular_weight()
    ip=X.isoelectric_point()
    extinction_coefficient=X.molar_extinction_coefficient()[0]

    ## Check the number of W to determine if the extinction coefficient is faithful

    return mw, ip, extinction_coefficient 


def check_enzyme_cut(enzyme, sequence):
    '''
    This function checks possible enzyme restriction cutting sites and replaces problematic regions with equivalent sequences.

    Input:
    enzyme --> List of enzymes to check.
    sequence --> DNA sequence to check.

    Output:
    new_sequence --> Sequence without the restriction sites.
    '''
    amino_acid_codons = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],  # Alanine
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # Arginine
        'N': ['AAT', 'AAC'],  # Asparagine
        'D': ['GAT', 'GAC'],  # Aspartic acid
        'C': ['TGT', 'TGC'],  # Cysteine
        'Q': ['CAA', 'CAG'],  # Glutamine
        'E': ['GAA', 'GAG'],  # Glutamic acid
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],  # Glycine
        'H': ['CAT', 'CAC'],  # Histidine
        'I': ['ATT', 'ATC', 'ATA'],  # Isoleucine
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],  # Leucine
        'K': ['AAA', 'AAG'],  # Lysine
        'M': ['ATG'],  # Methionine (Start codon)
        'F': ['TTT', 'TTC'],  # Phenylalanine
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],  # Proline
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # Serine
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],  # Threonine
        'W': ['TGG'],  # Tryptophan
        'Y': ['TAT', 'TAC'],  # Tyrosine
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],  # Valine
        '*': ['TAA', 'TAG', 'TGA']  # Stop codons
    }
    restriction_enzymes = {
        'EcoRI': 'GAATTC',
        'BamHI': 'GGATCC',
        'HindIII': 'AAGCTT',
        'NotI': 'GCGGCCGC',
        'XhoI': 'CTCGAG',
        'PstI': 'CTGCAG',
        'SacI': 'GAGCTC',
        'KpnI': 'GGTACC',
        'SmaI': 'CCCGGG',
        'XbaI': 'TCTAGA',
        'SpeI': 'ACTAGT',
        'NcoI': 'CCATGG',
        'SalI': 'GTCGAC',
        'ApaI': 'GGGCCC',
        'HaeIII': 'GGCC',
        'AluI': 'AGCT',
        'TaqI': 'TCGA',
        'BglII': 'AGATCT',
        'ClaI': 'ATCGAT',
        'MluI': 'ACGCGT',
        'BsaI': 'GGTCTC'
    }
    # Get restriction sequences for selected enzymes
    selected_re = [restriction_enzymes[y] for y in enzyme]
    # Convert sequence to a mutable list
    list_sequence = list(sequence.upper())
    # Iterate through the sequence
    for i in range(len(sequence) - 5):  # Ensuring enough length for restriction sites
        sliding_window = sequence[i:i+6]  # Most restriction sites are 6 bp long
        found_enzyme = next((enzyme_re for enzyme_re in selected_re if enzyme_re in sliding_window), None)
        if found_enzyme is not None:
            problematic_enzyme=list(restriction_enzymes.keys())[list(restriction_enzymes.values()).index(found_enzyme)]
            print(f'Sequence of cut for {problematic_enzyme} has been found')
            #Find the codon at position i+6
            codon_index=math.floor((i+3)/3) #zero indexed
            codon_identity=sequence[codon_index*3:codon_index*3+3]   
            # Find the amino acid corresponding to the codon
            residue = None
            for aa, codons in amino_acid_codons.items():
                if codon_identity in codons:
                    residue = aa
                    if residue in ['M', 'W']:
                        codon_index-=1
                        codon_identity=sequence[codon_index*3:codon_index*3+3]
                        for aa, codons in amino_acid_codons.items():
                            if codon_identity in codons:
                                residue = aa
                    print(f'Searching for a new codon for {residue}')
                    break
            if residue != '*':
                # Select an alternative codon (excluding the original one)
                alternative_codons = [codon for codon in amino_acid_codons[residue] if codon != codon_identity] 
                if alternative_codons:
                    alternative_codon = random.choice(alternative_codons)  # Pick a random alternative codon
                    list_sequence[codon_index*3:codon_index*3+3] = alternative_codon  # Replace codon
                    print(f'Replacing {"".join(codon_identity)} with {"".join(alternative_codon)} ')
                    sequence=''.join(list_sequence)
                else:
                    print('M and W next to each other, regenerating DNA seq to avoid restriction site')
                    return '',False
            else:
                print('The codon selected is a stop codon, must not be modified, redesigning DNA seq') # Probably I should select other residue
                return '',False
    # Convert list back to string
    new_sequence = ''.join(list_sequence)
    return new_sequence.lower(),True

