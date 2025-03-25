import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np
from Bio import PDB
import pyrosetta
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover


def plot_distributions_and_correlations(csv_file, reference_values=None):
    df = pd.read_csv(csv_file)

    # Extract the label from path
    df["label"] = df["path"].str.split("_").str[0]

    # Define color palette
    unique_labels = df["label"].unique()
    palette = sns.color_palette("husl", len(unique_labels))  # Different colors for each label
    label_colors = dict(zip(unique_labels, palette))

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))  # 2 rows, 3 columns
    metrics = ['dG', 'dSASA', 'sc']

    # Plot distributions (histograms)
    for i, metric in enumerate(metrics):
        for label in unique_labels:
            subset = df[df["label"] == label]
            sns.histplot(subset[metric], kde=True, ax=axes[0, i], label=label, color=label_colors[label], alpha=0.6)
        
        axes[0, i].set_title(f'Distribution of {metric}')
        axes[0, i].set_xlabel(metric)
        axes[0, i].set_ylabel('Count')
        axes[0, i].legend(title="Label")

        if reference_values and metric in reference_values:
            axes[0, i].axvline(x=reference_values[metric], color='black', linestyle='--', linewidth=2, label='Reference')

    # Plot correlations (scatter plots)
    scatter_plots = [('sc', 'dG'), ('dSASA', 'dG'), ('sc', 'dSASA')]
    for i, (x_metric, y_metric) in enumerate(scatter_plots):
        for label in unique_labels:
            subset = df[df["label"] == label]
            sns.scatterplot(x=x_metric, y=y_metric, data=subset, ax=axes[1, i], label=label, color=label_colors[label], alpha=0.6)

        axes[1, i].set_title(f'{y_metric} vs {x_metric}')
        axes[1, i].set_xlabel(x_metric)
        axes[1, i].set_ylabel(y_metric)

        # Fit and plot regression line
        x, y = df[x_metric], df[y_metric]
        m, b = np.polyfit(x, y, 1)
        axes[1, i].plot(x, m*x + b, color='red', linestyle='--')

        # Compute and display Pearson correlation (r) and R²
        r = df[y_metric].corr(df[x_metric])
        r2 = r**2
        axes[1, i].annotate(f'r = {r:.3f}\nR² = {r2:.3f}', xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10,
                            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray', alpha=0.8))

        if reference_values and x_metric in reference_values and y_metric in reference_values:
            axes[1, i].axvline(x=reference_values[x_metric], color='black', linestyle='--', linewidth=1)
            axes[1, i].axhline(y=reference_values[y_metric], color='black', linestyle='--', linewidth=1)

    plt.tight_layout()
    plt.savefig('combined_figure.png')
    plt.show()

    print("Combined figure with distributions and correlations saved as 'combined_figure.png'")


def substitute_residues_with_glycine(pdb_file, output_pdb, chain_id):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    io = PDB.PDBIO()
    
    class GlycineMutator(PDB.Select):
        def accept_residue(self, residue):
            if residue.parent.id == chain_id:
                # Change the residue name to GLY
                residue.resname = "GLY"
                
                # Remove all side chain atoms (keep backbone: N, CA, C, O)
                keep_atoms = {"N", "CA", "C", "O"}
                atoms_to_remove = [atom for atom in residue if atom.name not in keep_atoms]
                for atom in atoms_to_remove:
                    residue.detach_child(atom.id)
                    
            return True
    
    io.set_structure(structure)
    io.save(output_pdb, GlycineMutator())

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
        return result

    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot distributions of metrics with optional reference values')
    parser.add_argument('--csv', type=str, default='interface_metrics.csv', help='Path to the CSV file (default: interface_metrics.csv)')
    parser.add_argument('--ref', required=True, type=str, help='Reference PDB to compare with')
    args = parser.parse_args()
    
    ref_pdb = substitute_residues_with_glycine(args.ref, 'temp.pdb', 'A')
    ref_metrics = compute_interface_metrics('temp.pdb')
    
    plot_distributions_and_correlations(args.csv, ref_metrics)
