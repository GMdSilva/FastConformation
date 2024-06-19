import json
import os
import shutil

import pandas as pd

from glob import glob

import MDAnalysis as mda
from MDAnalysis.analysis import align


def parabola(x, a, b, c):
    return a * x**2 + b * x + c


def create_directory(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


def load_config(config_file):
    try:
        with open(config_file, 'r') as file:
            config = json.load(file)
        return config
    except FileNotFoundError:
        print(f"File {config_file} not found.")
    except json.JSONDecodeError:
        print(f"Error decoding JSON from the file {config_file}.")


def load_frames(file_list):
    """ Load frames from a list of files """
    universes = [mda.Universe(f) for f in file_list]
    return universes


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
    print(f"Loaded {len(u.trajectory)} predictions from {folder_path}")
    return u


def load_predictions(predictions_path, seq_pairs, jobname):
    predictions_dict = {}
    for seq_pair in seq_pairs:
        max_seq = seq_pair[0]
        extra_seq = seq_pair[1]
        folder_path = f"{predictions_path}/{jobname}_{max_seq}_{extra_seq}"

        universe = load_pdb_files_as_universe(folder_path)

        params = {'max_seq': max_seq,
                  'extra_seq': extra_seq,
                  'jobname': jobname,
                  'mda_universe': universe}
        predictions_dict[f'{jobname}_{max_seq}_{extra_seq}'] = params
    return predictions_dict


def load_predictions_json(predictions_path, seq_pairs, jobname):
    plddt_dict = {}

    for seq_pair in seq_pairs:
        max_seq = seq_pair[0]
        extra_seq = seq_pair[1]
        folder_path = f"{predictions_path}/{jobname}_{max_seq}_{extra_seq}"
        json_files = glob(f"{folder_path}/*.json")
        all_plddts = []

        for json_file in json_files:
            plddt_data = load_config(json_file)
            per_residue_plddt = plddt_data.get('plddt')
            if per_residue_plddt:
                all_plddts.append(per_residue_plddt)

        partial_plddt_dict = {'max_seq': max_seq,
                              'extra_seq': extra_seq,
                              'jobname': jobname,
                              'all_plddts': all_plddts}

        plddt_dict[f'{jobname}_{max_seq}_{extra_seq}'] = partial_plddt_dict

    return plddt_dict


def reorder_frames_by(frames, values):
    top = frames._topology
    frame_value_pairs = list(enumerate(values))
    sorted_values = sorted(frame_value_pairs, key=lambda x: x[1])  # Sort by rmsd_value (x[1])
    sorted_indices = [idx for idx, _ in sorted_values]
    traj_temp = []
    create_directory('tmp_frames')

    for idx, frame in enumerate(frames.trajectory[sorted_indices]):
        frame_name = f"tmp_frames/tmp_{idx}.pdb"
        ag = frames.select_atoms('all')
        ag.write(frame_name)
        traj_temp.append(frame_name)

    traj = mda.Universe(top, traj_temp, dt=1)
    align.AlignTraj(traj, traj, select='backbone', ref_Frame=0, in_memory=True).run()

    return traj


def save_traj(universe, traj_output_path, jobname, max_seq, extra_seq, traj_format, ordered):
    ag = universe.select_atoms('all')

    if not ordered:
        ordered = 'raw'

    ag.write(f"{traj_output_path}/{jobname}_{max_seq}_{extra_seq}_ordered_by_{ordered}.{traj_format}", frames='all')
    shutil.rmtree('tmp_frames')


def auto_select_2d_references(references_dataset_path, analysis_type):
    df_references = pd.read_csv(references_dataset_path)
    # Get unique values in the 'trial' column
    unique_trials = df_references['trial'].unique()
    all_modes_diff = {}
    modes_diff_partial = {}
    # Iterate over each unique trial and create a subset
    for trial in unique_trials:
        subset = df_references[df_references['trial'] == trial]
        modes_subset = list(subset[analysis_type])
        if len(modes_subset) > 1:
            max_difference = 0
            index1 = 0
            index2 = 0
            for i in range(len(modes_subset)):
                for j in range(i + 1, len(modes_subset)):
                    difference = abs(modes_subset[i] - modes_subset[j])
                    if difference > max_difference:
                        max_difference = difference
                        index1 = i
                        index2 = j
                modes_diff_partial = {"trial": trial,
                                      "max_difference": max_difference,
                                      "index1": index1,
                                      "index2": index2}
        else:
            modes_diff_partial = {"trial": trial,
                                  "max_difference": None,
                                  "index1": None,
                                  "index2": None}
        all_modes_diff[trial] = modes_diff_partial

    max_diff_key = None
    max_diff_value = float('-inf')

    # Iterate over the dictionary items
    for key, value in all_modes_diff.items():
        # Skip entries with None as max_difference
        if value['max_difference'] is not None and value['max_difference'] > max_diff_value:
            max_diff_value = value['max_difference']
            max_diff_key = key

    ref1_path = df_references.loc[df_references['trial'] == max_diff_key, 'pdb_filename'].values[all_modes_diff[max_diff_key]['index1']]
    ref2_path = df_references.loc[df_references['trial'] == max_diff_key, 'pdb_filename'].values[all_modes_diff[max_diff_key]['index2']]

    return ref1_path, ref2_path