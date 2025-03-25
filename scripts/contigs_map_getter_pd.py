#!/usr/bin/env python3
'''
Script to calculate the contigs.contigsmap that RFD needs.

Input:

-i|--input --> pdb file of the input

Output:

Prints in the screen the groupped residue number in the contigsmap format

'''
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', help='protein pdbfile path', required=False)
args = parser.parse_args()



def parse_pdb(pdb_file):
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    chain_A_residues = []
    chain_B_residues = []

    for line in lines:
        if line.startswith('ATOM'):
            chain = line[21]
            residue_id = int(line[22:26].strip())
            if chain == 'A':
                chain_A_residues.append(residue_id)
            elif chain == 'B':
                chain_B_residues.append(residue_id)

    residue_length_chain_A = len(set(chain_A_residues))
    residue_length_chain_B = len(set(chain_B_residues))

    residue_start_B = min(chain_B_residues)
    residue_finish_B = max(chain_B_residues)

    return f"[ {residue_length_chain_A}-{residue_length_chain_A}/0 B{residue_start_B}-{residue_finish_B} ]"

if __name__ == "__main__":
    pdb_file = args.input
    result = parse_pdb(pdb_file)
    print(result)
