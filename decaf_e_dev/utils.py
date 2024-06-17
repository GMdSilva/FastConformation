import json
import os
import shutil


from glob import glob

import MDAnalysis as mda


def parabola(x, a, b, c):
    return a * x**2 + b * x + c


def create_directory(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)
    print(f"Directory '{path}' created successfully.")


def load_config(config_file):
    try:
        with open(config_file, 'r') as file:
            config = json.load(file)
        print(f"Successfully loaded config file {config_file}.")
        return config
    except FileNotFoundError:
        print(f"Configuration file {config_file} not found.")
    except json.JSONDecodeError:
        print(f"Error decoding JSON from the configuration file {config_file}.")


def load_frames(file_list):
    """ Load frames from a list of files """
    universes = [mda.Universe(f) for f in file_list]
    return universes


def save_reordered_trajectory(ordered_file_list, output_file):
    """ Save the reordered frames as a new trajectory """
    with mda.Writer(output_file, n_atoms=universes[0].atoms.n_atoms) as W:
        for pdb in ordered_file_list:
            u = mda.Universe(pdb)
            W.write(u.atoms)


def load_pdb_files_as_universe(folder_path):
    """
    Load all PDB files in the specified folder as a Universe,
    using the first PDB file (sorted alphabetically) as the topology.

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
    u = mda.Universe(topology, *pdb_files, dt=1)

    # Print some information about the loaded universe
    print(f"Loaded Universe with {len(u.trajectory)} frames at {folder_path}.")
    return u


def load_predictions(output_path, seq_pairs, jobname):
    predictions_dict = {}
    for seq_pair in seq_pairs:
        max_seq = seq_pair[0]
        extra_seq = seq_pair[1]
        folder_path = f"{output_path}/{jobname}_{max_seq}_{extra_seq}"

        universe = load_pdb_files_as_universe(folder_path)

        params = {'max_seq': max_seq,
                  'extra_seq': extra_seq,
                  'jobname': jobname,
                  'mda_universe': universe}
        predictions_dict[f'{jobname}_{max_seq}_{extra_seq}'] = params
    return predictions_dict



def load_predictions_json(output_path, seq_pairs, jobname):
    plddt_dict = {}
    for seq_pair in seq_pairs:
        max_seq = seq_pair[0]
        extra_seq = seq_pair[1]
        folder_path = f"{output_path}/{jobname}_{max_seq}_{extra_seq}"

        json_files = glob(f"{folder_path}/*.json")

        print(f"Reading {len(json_files)} json metadata files at {folder_path}.")

        all_plddts = {}
        for json_file in json_files:
            plddt_data = load_config(json_file)
            per_residue_plddt = plddt_data.get('plddt')
            if per_residue_plddt:
                all_plddts[json_file] = per_residue_plddt

        partial_plddt_dict = {'max_seq': max_seq,
                      'extra_seq': extra_seq,
                      'jobname': jobname,
                      'all_plddts': all_plddts}

        plddt_dict[f'{jobname}_{max_seq}_{extra_seq}'] = partial_plddt_dict
    return plddt_dict



