#! /usr/bin/env python3
"""
Script to redesign pdb for Partial Diffusion using BioPython
"""

import argparse
from Bio import PDB

def modify_pdb(input_pdb, output_pdb):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    new_structure = []
    prev_resid = None

    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                resid = residue.get_id()[1]

                if chain_id == 'A':
                    prev_resid = resid
                    new_structure.append(residue)
                elif chain_id == 'B':
                    if prev_resid is not None and resid != prev_resid + 1:
                        residue.id = (' ', prev_resid + 1, ' ')
                        prev_resid += 1
                    new_structure.append(residue)

    # Writing the modified structure to a new PDB file
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_file", required=True, help="Input PDB file")
    args = parser.parse_args()

    output_pdb = f"{args.pdb_file.split('.')[0]}_mod.pdb"
    modify_pdb(args.pdb_file, output_pdb)

    print(f"Modified PDB saved as {output_pdb}")
