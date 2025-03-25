#!/usr/bin/env python3

import os
import subprocess
import argparse

# Functions

def pdb_extract_seq(input_path, chain):
    """Extract sequences from PDB files."""
    command = f'python3 /apps/scripts/protein_design/scripts/pdb_extract_seq.py -i {input_path} -o {input_path}/hits.fasta -c {chain}'
    subprocess.run(command, shell=True)

def add_initial_met(fasta_path):
    """Add an initial Met (M) if not already present."""
    input_fasta = os.path.join(fasta_path, "hits.fasta")
    output_fasta = os.path.join(fasta_path, "hits_met.fasta")

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                outfile.write(line)  # Write the header line as is
            else:
                line = line.strip()
                if not line.startswith('M'):
                    outfile.write('M' + line + 'SGG\n')  # Add 'M' and 'SGG'
                else:
                    outfile.write(line + 'SGG\n')  # Just add 'SGG'


def convert_fasta_to_csv(input_fasta, output_csv):
    """Convert the final fasta file to a CSV format."""
    with open(input_fasta, 'r') as fasta_file, open(output_csv, 'w') as csv_file:
        name, sequence = "", ""
        for line in fasta_file:
            if line.startswith('>'):
                if sequence:
                    csv_file.write(f"{name},{sequence}\n")
                name = line[1:].strip()  # Remove the ">" and strip spaces
                sequence = ""
            else:
                sequence += line.strip()
        if sequence:  # Write the last sequence
            csv_file.write(f"{name},{sequence}\n")

# Main

def main():
    parser = argparse.ArgumentParser(description="Automate sequence extraction and optimization")
    parser.add_argument("-i", "--input", required=True, help="Path to input files")
    parser.add_argument("-c", "--chain", required=True, help="Chain identifier for sequence extraction")
    args = parser.parse_args()

    input_path = args.input
    chain = args.chain

    # Step 1: Extract sequences
    pdb_extract_seq(input_path, chain)
    
    # Step 2: Add an initial Met if necessary and final GGS
    add_initial_met(input_path)
    
    # Step 3: Convert to CSV
    convert_fasta_to_csv(f"{input_path}/hits_met.fasta", f"{input_path}/hits_met.csv")

#    overhang_5 = "ttttgtttaactttaataaggagatatacc" # EGFP-10xHis(noATG) /pRSFDuet_Peptide1-His
#    overhang_3 = "ttggaggtcttgtttcagggt" # EGFP-10xHis(noATG) / pRSFDuet_Peptide1-His



if __name__ == "__main__":
    main()
