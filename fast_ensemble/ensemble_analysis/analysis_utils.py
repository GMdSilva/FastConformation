import json
import os
import shutil
import subprocess
import pandas as pd
from glob import glob
import MDAnalysis as mda
from MDAnalysis.analysis import align

def parabola(x, a, b, c):
    """
    Calculate the value of a parabola given the coefficients.

    Parameters:
    x (float or array-like): The independent variable.
    a (float): Coefficient for the quadratic term.
    b (float): Coefficient for the linear term.
    c (float): Constant term.

    Returns:
    float or array-like: The value of the parabola at x.
    """
    return a * x**2 + b * x + c


def create_directory(path):
    """
    Create a directory at the specified path. If the directory already exists, it is removed and recreated.

    Parameters:
    path (str): The path to the directory.

    Raises:
    OSError: If the directory cannot be created.
    """
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


def load_config(config_file):
    """
    Load a JSON configuration file.

    Parameters:
    config_file (str): The path to the configuration file.

    Returns:
    dict: The configuration as a dictionary.

    Raises:
    FileNotFoundError: If the file does not exist.
    json.JSONDecodeError: If the file is not a valid JSON.
    """
    try:
        with open(config_file, 'r') as file:
            config = json.load(file)
        return config
    except FileNotFoundError:
        print(f"File {config_file} not found.")
    except json.JSONDecodeError:
        print(f"Error decoding JSON from the file {config_file}.")


def load_frames(file_list):
    """
    Load a list of molecular dynamics files as MDAnalysis Universes.

    Parameters:
    file_list (list): A list of file paths to load.

    Returns:
    list: A list of MDAnalysis Universes.
    """
    universes = [mda.Universe(f) for f in file_list]
    return universes


def load_pdb_files_as_universe(folder_path, reindex):
    """
    Load all PDB files in the specified folder as a Universe, using the first PDB file as the topology.

    Parameters:
    folder_path (str): Path to the folder containing PDB files.
    reindex (int or None): If provided, reindex the PDB files so the first residue matches this index.

    Returns:
    MDAnalysis.Universe: The loaded Universe.

    Raises:
    FileNotFoundError: If no PDB files are found in the folder.
    """
    pdb_files = sorted(glob(os.path.join(folder_path, '*.pdb')))

    if not pdb_files:
        raise FileNotFoundError("No PDB files found in the specified folder")

    topology = pdb_files[0]

    if reindex:
        temp_pdb_path = os.path.join(folder_path, 'temp.pdb')
        command = f"pdb_reres -{reindex} \"{topology}\" > \"{temp_pdb_path}\""
        subprocess.run(command, shell=True, check=True)
        topology = temp_pdb_path

    u = mda.Universe(topology, *pdb_files, dt=1)

    if reindex:
        os.remove(temp_pdb_path)

    print(f"Loaded {len(u.trajectory)} predictions from {folder_path}")
    return u


def load_predictions(predictions_path, seq_pairs, jobname, starting_residue):
    """
    Load predictions from a set of PDB files as MDAnalysis Universes.

    Parameters:
    predictions_path (str): Path to the directory containing the predictions.
    seq_pairs (list): A list of sequence pairs for the predictions.
    jobname (str): The job name associated with the predictions.
    starting_residue (int or None): The starting residue index for reindexing.

    Returns:
    dict: A dictionary of Universes and associated metadata.
    """
    predictions_dict = {}
    
    for seq_pair in seq_pairs:
        max_seq, extra_seq = seq_pair
        folder_pattern = os.path.join(predictions_path, f"*{max_seq}_{extra_seq}*")

        # Find all folders that match the pattern
        matching_folders = [f for f in glob(folder_pattern) if os.path.isdir(f)]
        
        # Check if any matching folders are found
        if not matching_folders:
            print(f"Warning: No folders found for sequence pair '{max_seq}_{extra_seq}'. Skipping...")
            continue

        for folder_path in matching_folders:
            # Search for files containing max_seq_extra_seq within the matching folder
            file_pattern = os.path.join(folder_path, f"*{max_seq}_{extra_seq}*.pdb")
            pdb_files = glob(file_pattern)

            # Check if any matching files are found
            if not pdb_files:
                print(f"Warning: No files containing '{max_seq}_{extra_seq}' found in folder '{folder_path}'. Skipping...")
                continue

            # Load the found PDB files as MDAnalysis Universes
            universe = load_pdb_files_as_universe(pdb_files, starting_residue)

            params = {'max_seq': max_seq,
                      'extra_seq': extra_seq,
                      'jobname': jobname,
                      'mda_universe': universe}

            predictions_dict[f'{jobname}_{max_seq}_{extra_seq}'] = params
    
    return predictions_dict

def load_predictions_json(predictions_path, seq_pairs, jobname):
    """
    Load pLDDT scores from JSON files associated with predictions.

    Parameters:
    predictions_path (str): Path to the directory containing the predictions.
    seq_pairs (list): A list of sequence pairs for the predictions.
    jobname (str): The job name associated with the predictions.

    Returns:
    dict: A dictionary of pLDDT scores for each prediction.
    """
    plddt_dict = {}

    for seq_pair in seq_pairs:
        max_seq, extra_seq = seq_pair
        folder_pattern = os.path.join(predictions_path, f"*{max_seq}_{extra_seq}*")

        # Find all folders that match the pattern
        matching_folders = [f for f in glob(folder_pattern) if os.path.isdir(f)]

        # Check if any matching folders are found
        if not matching_folders:
            print(f"Warning: No folders found for sequence pair '{max_seq}_{extra_seq}'. Skipping...")
            continue

        all_plddts = []

        for folder_path in matching_folders:
            # Search for JSON files containing max_seq_extra_seq within the matching folder
            json_file_pattern = os.path.join(folder_path, f"*{max_seq}_{extra_seq}*.json")
            json_files = glob(json_file_pattern)

            # Check if any matching files are found
            if not json_files:
                print(f"Warning: No JSON files containing '{max_seq}_{extra_seq}' found in folder '{folder_path}'. Skipping...")
                continue

            for json_file in json_files:
                # Load pLDDT data from the JSON file
                plddt_data = load_config(json_file)
                per_residue_plddt = plddt_data.get('plddt')

                if per_residue_plddt:
                    all_plddts.append(per_residue_plddt)

        if all_plddts:  # Only add to dictionary if there are valid pLDDT values
            partial_plddt_dict = {
                'max_seq': max_seq,
                'extra_seq': extra_seq,
                'jobname': jobname,
                'all_plddts': all_plddts
            }

            plddt_dict[f'{jobname}_{max_seq}_{extra_seq}'] = partial_plddt_dict
        else:
            print(f"Warning: No valid pLDDT data found for sequence pair '{max_seq}_{extra_seq}'. Skipping...")

    return plddt_dict


def reorder_frames_by(frames, values):
    """
    Reorder trajectory frames based on associated values (e.g., RMSD).

    Parameters:
    frames (MDAnalysis.Universe): The Universe containing the trajectory frames to reorder.
    values (list): A list of values associated with each frame.

    Returns:
    MDAnalysis.Universe: A new Universe with frames reordered by the associated values.
    """
    top = frames._topology
    frame_value_pairs = list(enumerate(values))
    sorted_values = sorted(frame_value_pairs, key=lambda x: x[1])
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
    """
    Save the trajectory of a Universe to a file.

    Parameters:
    universe (MDAnalysis.Universe): The Universe containing the trajectory to save.
    traj_output_path (str): The directory where the trajectory file will be saved.
    jobname (str): The job name associated with the trajectory.
    max_seq (str): The maximum sequence associated with the trajectory.
    extra_seq (str): The extra sequence associated with the trajectory.
    traj_format (str): The format to save the trajectory in (e.g., 'pdb').
    ordered (str): A description of the ordering of the frames (default is 'raw').

    Returns:
    None
    """
    ag = universe.select_atoms('all')

    if not ordered:
        ordered = 'raw'

    ag.write(f"{traj_output_path}/{jobname}_{max_seq}_{extra_seq}_ordered_by_{ordered}.{traj_format}", frames='all')
    shutil.rmtree('tmp_frames')


def auto_select_2d_references(references_dataset_path, analysis_type):
    """
    Automatically select two 2D references based on the most distant mode representatives.

    Parameters:
    references_dataset_path (str): Path to the CSV file containing reference data.
    analysis_type (str): The type of analysis (e.g., 'tmscore', 'rmsd') to use for mode selection.

    Returns:
    tuple: Paths to the two selected reference structures.
    """
    df_references = pd.read_csv(references_dataset_path)
    unique_trials = df_references['trial'].unique()
    all_modes_diff = {}

    for trial in unique_trials:
        subset = df_references[df_references['trial'] == trial]
        modes_subset = list(subset[analysis_type])
        modes_diff_partial = {}

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

    for key, value in all_modes_diff.items():
        if value['max_difference'] is not None and value['max_difference'] > max_diff_value:
            max_diff_value = value['max_difference']
            max_diff_key = key

    ref1_path = df_references.loc[df_references['trial'] == max_diff_key, 'pdb_filename'].values[all_modes_diff[max_diff_key]['index1']]
    ref2_path = df_references.loc[df_references['trial'] == max_diff_key, 'pdb_filename'].values[all_modes_diff[max_diff_key]['index2']]

    return ref1_path, ref2_path
