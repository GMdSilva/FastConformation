import os
import MDAnalysis as mda
from glob import glob


def load_pdb_files_as_universe(folder_path):
    """
    Load all PDB files in the specified folder as a Universe, using the first PDB file (sorted alphabetically) as the topology.

    Parameters:
    folder_path (str): Path to the folder containing PDB files.

    Returns:
    MDAnalysis.Universe: The loaded Universe.
    """
    # Get a list of all PDB files in the folder, sorted alphabetically
    pdb_files = sorted(glob(os.path.join(folder_path, '*.pdb')))

    # Check if there are any PDB files in the folder
    if not pdb_files:
        raise FileNotFoundError("No PDB files found in the specified folder.")

    # Use the first PDB file as the topology
    topology = pdb_files[0]

    # Load all PDB files as a Universe using the first file as the topology
    u = mda.Universe(topology, *pdb_files)

    # Print some information about the loaded universe
    print(f"Loaded Universe with {len(u.trajectory)} frames")
    print(f"Using topology from: {topology}")

    # Example: Accessing atoms and residues
    print(f"Number of atoms: {len(u.atoms)}")
    print(f"Number of residues: {len(u.residues)}")

    return u


# Example usage
folder_path = 'decaf_e_dev/data/sample_predictions/abl1/abl_wt_64_128'
universe = load_pdb_files_as_universe(folder_path)
print(universe.atoms)
