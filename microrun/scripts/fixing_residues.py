#!/usr/bin/env python3 
'''
Script to fix residues in the case for sequence diversity. Uses the function define in biopython_align.py

Input: 

--fixed: List of residues to be fixed. Formatted as "[1,3,10-20]
--pdb_input: Path to the pdb input

Output:

The PDB inputted is modified so after the last TER remarks fixing each of the defined residues is added
'''


import argparse
from biopython_align import add_fixed_residues

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument("--residues", required=False, default=None, help='List of residues of the binder to fix (Useful for scaffolding)')
    parser.add_argument("--pdb_input", required=True, help='Path to the pdb file')
    args=parser.parse_args()

    add_fixed_residues(output_path=args.pdb_input, residues=args.residues)
