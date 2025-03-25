#!/usr/bin/env python3
'''
Compute relevant Rosetta statistics for RFD backbones: shape complementarity, dG and dSASA.
Input: A list of PDBs to run on.
Output: interface_metrics.csv, containing parameters and path.
Usage: python scoring_backbones.py "output*/design_*pdb" [--processes N]
'''
import pyrosetta
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
import pandas as pd
import glob
import sys
import multiprocessing
import argparse
import os
from functools import partial

def compute_interface_metrics(pdb_path, chain1='A', chain2='B'):
    # Initialize pyrosetta in each worker process
    # Setting quiet=True to reduce output spam in parallel execution
    pyrosetta.init(options="-mute all")
    
    try:
        pose = pyrosetta.pose_from_pdb(pdb_path)
        interface = f'{chain1}_{chain2}'
        ia_mover = InterfaceAnalyzerMover(interface)
        ia_mover.apply(pose)
        interface_dG = ia_mover.get_interface_dG()
        interface_dSASA = ia_mover.get_interface_delta_sasa()
        sc = ia_mover.get_all_data().sc_value
        
        result = {
            "path": pdb_path,
            "dG": round(interface_dG, 1),
            "dSASA": round(interface_dSASA, 1),
            "sc": round(sc, 3)
        }
        print(f"Processed {os.path.basename(pdb_path)}")
        return result
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return None

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Compute Rosetta interface metrics in parallel')
    parser.add_argument('file_pattern', help='Glob pattern for PDB files (e.g., "output*/design_*pdb")')
    parser.add_argument('--processes', '-p', type=int, default=4,
                        help=f'Number of parallel processes (default: 4)')
    parser.add_argument('--chains', '-c', default='A_B', help='Interface chains (format: A_B)')
    args = parser.parse_args()
    
    # Parse chains
    chains = args.chains.split('_')
    if len(chains) != 2:
        print("Error: Chain format should be 'A_B'")
        sys.exit(1)
    chain1, chain2 = chains
    
    # Find PDB files
    print(f"Searching for files matching {args.file_pattern}")
    pdb_files = glob.glob(args.file_pattern)
    
    if not pdb_files:
        print("No matching files found.")
        sys.exit(1)
    else:
        print(f"Found {len(pdb_files)} files.")
        print(f"Using {args.processes} processes for parallel execution.")
    
    # Create a pool of worker processes
    pool = multiprocessing.Pool(processes=args.processes)
    
    # Create a partial function with the chain arguments
    worker_func = partial(compute_interface_metrics, chain1=chain1, chain2=chain2)
    
    # Process files in parallel
    results = pool.map(worker_func, pdb_files)
    
    # Close the pool and wait for work to finish
    pool.close()
    pool.join()
    
    # Filter out None results (failed processing)
    results = [r for r in results if r is not None]
    
    # Create dataframe and save results
    df = pd.DataFrame(results)
    df.to_csv("interface_metrics.csv", index=False)
    
    print(f"Results saved to interface_metrics.csv ({len(results)} of {len(pdb_files)} files processed successfully)")

if __name__ == "__main__":
    main()