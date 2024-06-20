import subprocess
import os
import re
import shutil

import pandas as pd
import numpy as np
from glob import glob

import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

import MDAnalysis as mda

from tqdm import tqdm

from decaf_e_dev.ensemble_analysis.analysis_utils import create_directory

TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'


def slice_models(universe, selection, temp_traj_path):
    for idx, ts in enumerate(universe.trajectory):
        ag = universe.select_atoms(selection)
        ag.write(f"{temp_traj_path}/tmp_{idx:0>5}.pdb")


def tmscore_wrapper(mobile, reference):
    command = f"TMscore {mobile} {reference} -outfmt 2"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    tmscore = None

    # Read line by line as they are output
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            match = re.search(r'TM-score\s+=\s+([0-9.]+)', output)
            if match:
                tmscore = float(match.group(1))
    process.poll()
    return tmscore


def run_tmscore(folder_path, custom_ref):
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
    if not custom_ref:
        custom_ref = pdb_files[0]
    tmscore_dict = {}
    tmscores = []
    frames = []
    for idx, pdb_file in enumerate(pdb_files):
        frames.append(idx)
        tmscores.append(tmscore_wrapper(pdb_file, custom_ref))
    tmscore_dict['frame'] = frames
    tmscore_dict['tmscore'] = tmscores
    return tmscore_dict


def tmscore_kde(tmscore_data: list, input_dict: dict, slice_predictions) -> dict:
    jobname = input_dict['jobname']
    max_seq = input_dict['max_seq']
    extra_seq = input_dict['extra_seq']
    output_path = input_dict['output_path']
    # Calculate KDE
    kde = gaussian_kde(tmscore_data, bw_method='silverman')
    x_vals = np.linspace(min(tmscore_data), max(tmscore_data), 1000)
    kde_vals = kde(x_vals)

    # Find modes (peaks)
    peaks, _ = find_peaks(kde_vals, prominence=1)
    modes = x_vals[peaks]

    # Find the two most distant modes
    if len(modes) > 1:
        distances = np.abs(np.subtract.outer(modes, modes))
        np.fill_diagonal(distances, 0)
        max_dist_indices = np.unravel_index(np.argmax(distances), distances.shape)
        most_distant_modes = modes[list(max_dist_indices)]
    else:
        most_distant_modes = modes

    # Mark the most distant modes
    modes_dict = {}
    mode_n = 0
    for mode_density, mode in zip(kde_vals[peaks], modes):
        array = np.asarray(tmscore_data).flatten()  # Flatten the array to ensure it's 1-dimensional
        target_value = np.asarray(mode).item()  # Ensure target_value is a scalar

        mode_index = (np.abs(array - target_value)).argmin()
        mode_results = {'mode_index': mode_index,
                        'mode_value': mode,
                        'mode_density': mode_density}

        mode_n += 1

        modes_dict[f'mode_{mode_n}'] = mode_results

    # Plot KDE and mark modes
    plt.figure(figsize=(6, 3))
    plt.plot(x_vals, kde_vals, label='KDE')
    plt.scatter(modes, kde_vals[peaks], color='red', zorder=5, label='Modes')

    for mode in most_distant_modes:
        plt.axvline(mode, color='blue', linestyle='--', label='Most Distant Mode')

    plt.title(f'{jobname} max_seq: {max_seq} extra_seq: {extra_seq}', fontsize=16)
    if slice_predictions:
        plt.suptitle(f'{slice_predictions}', fontsize=10)
    plt.xlabel('TM-Score', fontsize=14)
    plt.ylabel('Density', fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)  # Major ticks
    plt.tick_params(axis='both', which='minor', labelsize=12)  # Minor ticks (if any)
    plt.grid(False)

    # Show the plot
    plt.tight_layout()

    figure_full_path = f"{output_path}/{jobname}/analysis/tmscore_1d/{jobname}_{max_seq}_{extra_seq}_tmscore_kde.png"

    plt.savefig(figure_full_path, dpi=300)
    plt.close()

    return modes_dict


def tmscore_mode_analysis(prediction_dicts, input_dict, custom_ref, slice_predictions):
    jobname = input_dict['jobname']
    output_path = input_dict['output_path']
    predictions_path = input_dict['predictions_path']

    print('')
    with tqdm(total=len(prediction_dicts), bar_format=TQDM_BAR_FORMAT) as pbar:
        for prediction in prediction_dicts:
            pbar.set_description(f'Running 1D TM-Score analysis for {prediction}')
            universe = prediction_dicts[prediction]['mda_universe']
            max_seq = prediction_dicts[prediction]['max_seq']
            extra_seq = prediction_dicts[prediction]['extra_seq']
            input_dict['max_seq'] = max_seq
            input_dict['extra_seq'] = extra_seq

            if slice_predictions:
                create_directory(f'{predictions_path}/tmp_pdb')
                slice_models(universe, slice_predictions,  f'{predictions_path}/tmp_pdb')
                prediction_dicts[prediction]['tmscore_dict'] = run_tmscore(f'{predictions_path}/tmp_pdb', custom_ref)

            else:
                prediction_dicts[prediction]['tmscore_dict'] = run_tmscore(predictions_path, custom_ref)

            if slice_predictions:
                shutil.rmtree(f'{predictions_path}/tmp_pdb')


            tmscore_df = pd.DataFrame.from_dict(prediction_dicts[prediction]['tmscore_dict'], orient='columns')
            full_df_path = f"{output_path}/{jobname}/analysis/tmscore_1d/{jobname}_{max_seq}_{extra_seq}_tmscore_1d_df.csv"
            tmscore_df.to_csv(full_df_path)

            tmscore_modes = tmscore_kde(prediction_dicts[prediction]['tmscore_dict']['tmscore'],
                                        input_dict,
                                        slice_predictions)

            for mode in tmscore_modes:
                mode_index = tmscore_modes[mode]['mode_index']
                mode_value = int(tmscore_modes[mode]['mode_value']*100)
                frame_to_save = prediction_dicts[prediction]['tmscore_dict']['frame'][mode_index]
                universe.trajectory[int(frame_to_save)]
                # Save the frame to a file

                full_pdb_path = (f"{output_path}/"
                                 f"{jobname}/"
                                 f"analysis/"
                                 f"representative_structures/"
                                 f"tmscore_1d/"
                                 f"{jobname}_"
                                 f'{max_seq}_'
                                 f'{extra_seq}_'
                                 f'frame_'
                                 f'{mode_index}_'
                                 f'tmscore_'
                                 f'{mode_value}.pdb')


                with mda.Writer(full_pdb_path, universe.atoms.n_atoms) as W:
                    W.write(universe.atoms)

                tmscore_modes[mode]['pdb_filename'] = full_pdb_path

            prediction_dicts[prediction]['tmscore_dict']['tmscore_modes'] = tmscore_modes
            pbar.update(n=1)

    return prediction_dicts


def build_dataset_tmscore_modes(results_dict, input_dict):
    jobname = input_dict['jobname']
    output_path = input_dict['output_path']

    trials = []
    mode_labels = []
    mode_indexes = []
    mode_densities = []
    mode_values = []
    pdb_filenames = []

    for result in results_dict:
        mode_data = results_dict[result]['tmscore_dict']['tmscore_modes']
        for mode in mode_data:
            trials.append(result)
            mode_labels.append(mode)
            mode_indexes.append(results_dict[result]['tmscore_dict']['tmscore_modes'][mode]['mode_index'])
            mode_densities.append(results_dict[result]['tmscore_dict']['tmscore_modes'][mode]['mode_density'])
            mode_values.append(results_dict[result]['tmscore_dict']['tmscore_modes'][mode]['mode_value'])
            pdb_filenames.append(results_dict[result]['tmscore_dict']['tmscore_modes'][mode]['pdb_filename'])

    df = pd.DataFrame({'mode_label': mode_labels,
                       'representative_Frame': mode_indexes,
                       'tmscore': mode_values,
                       'mode_peak_density': mode_densities,
                       'trial': trials,
                       'pdb_filename': pdb_filenames
                       })

    full_df_path = f"{output_path}/{jobname}/analysis/tmscore_1d/{jobname}_tmscore_1d_analysis_results.csv"

    print(f"\nSaving {jobname} TM-Score 1D analysis results to {full_df_path}\n")

    print(f"Representative Structures for {jobname} TM-Score 1D analysis saved at {output_path}/{jobname}/analysis/"
          f"representative_structures/tmscore_1d/\n")

    df.to_csv(full_df_path, index=False)
