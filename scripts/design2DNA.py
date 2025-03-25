#!/usr/bin/env python3

import os
import subprocess
import argparse

# Functions

def pdb_extract_seq(input_path, chain):
    """Extract sequences from PDB files."""
    command = f'python3 pdb_extract_seq.py -i {input_path} -o {input_path}/output.fasta -c {chain}'
    subprocess.run(command, shell=True)

def add_initial_met(fasta_path):
    """Add an initial Met if not already present."""
    command = f"sed '/^[^>M]/ s/^/M/' {fasta_path}/output.fasta > {fasta_path}/output_met.fasta"
    subprocess.run(command, shell=True)

def codon_optimization(input_fasta, output_fasta, max_relax):
    """Perform codon optimization for E. coli expression."""
    command = f'codon_harmony --input {input_fasta} --output {output_fasta} --verbose 1 --host "Escherichia coli B" --max-relax {max_relax} --cycles 50 --inner-cycles 50 --restriction-enzymes'
    subprocess.run(command, shell=True)

def merge_optimized_sequences(input_path):
    """Merge all optimized sequences into one file."""
    command = f'cat {input_path}/output_met_RevTrans*fasta >> {input_path}/hits_RTseqs.fasta'
    subprocess.run(command, shell=True)

def add_cloning_overhangs(input_fasta, output_fasta):
    """Add cloning overhangs to the sequences."""
    overhang_5 = "ttttgtttaactttaataaggagatatacc" # EGFP-10xHis(noATG)
    overhang_3 = "ttggaggtcttgtttcagggt" # EGFP-10xHis(noATG)
    command = f"""awk '{{if (substr($0,1,1)==">") print "\\n"$0; else printf("%s", $0)}}' {input_fasta} | awk '/^>/ {{printf("%s\\n", $0); next;}} {{print "{overhang_5}"$0"{overhang_3}";}}' > {output_fasta}"""
    subprocess.run(command, shell=True)
    # Remove the first empty line
    subprocess.run(f"sed -i '1d' {output_fasta}", shell=True)

def convert_fasta_to_csv(input_fasta, output_csv):
    """Convert the final fasta file to a CSV format."""
    command = f"""awk '/^>/ {{if (seq) print name "," seq; name = substr($0, 2); seq = ""; next}} \
{{seq = seq $0}} END {{if (seq) print name "," seq}}' {input_fasta} > {output_csv}"""
    subprocess.run(command, shell=True)

# Main

def main():
    parser = argparse.ArgumentParser(description="Automate sequence extraction and optimization")
    parser.add_argument("-i", "--input", required=True, help="Path to input files")
    parser.add_argument("-c", "--chain", required=True, help="Chain identifier for sequence extraction")
    args = parser.parse_args()

    input_path = args.input
    chain = args.chain

    # Step 1: Extract sequences
    pdb_extract_seq(input_path, input_path, chain)
    
    # Step 2: Add an initial Met if necessary
    add_initial_met(input_path)
    
    # Step 3: Codon optimization
    max_relax_values = [0.15, 0.25, 0.3]
    for idx, max_relax in enumerate(max_relax_values):
        input_fasta = f"{input_path}/output_met.fasta" if idx == 0 else f"{input_path}/failed_seqs.fasta"
        output_fasta = f"{input_path}/output_met_RevTrans_{idx+1}.fasta"
        codon_optimization(input_fasta, output_fasta, max_relax)

    # Step 4: Merge all optimized sequences
    merge_optimized_sequences(input_path)
    
    # Step 5: Add cloning overhangs
    add_cloning_overhangs(f"{input_path}/hits_RTseqs.fasta", f"{input_path}/hits_RTseqs_overhang.fasta")
    
    # Step 6: Convert to CSV
    convert_fasta_to_csv(f"{input_path}/hits_RTseqs_overhang.fasta", f"{input_path}/hits_RTseqs_overhang.fasta.csv")

if __name__ == "__main__":
    main()
