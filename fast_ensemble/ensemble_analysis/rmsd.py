import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

import MDAnalysis as mda
from MDAnalysis.analysis import rms

from tqdm import tqdm
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore

TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'

def calculate_rmsd(u: mda.Universe,
                   ref: mda.Universe = None,
                   align_range: str = "backbone",
                   analysis_range: str = "backbone") -> dict:
    """
    Calculate the Root Mean Square Deviation (RMSD) of a given molecular dynamics trajectory.

    Parameters:
    u (mda.Universe): The MDAnalysis Universe containing the trajectory to be analyzed.
    ref (mda.Universe): The reference Universe for RMSD calculation. If None, the first frame of `u` is used.
    align_range (str): The atom selection string for alignment of the trajectory (default is "backbone").
    analysis_range (str): The atom selection string for RMSD calculation (default is "backbone").

    Returns:
    dict: A dictionary containing the RMSD values for each frame with keys 'frame', `align_range`, and `analysis_range`.
    """
    # if reference is not supplied, use first frame
    if ref is None:
        ref = u.select_atoms('protein')
    else:
        ref_u = mda.Universe(ref)
        ref = ref_u.select_atoms('protein')

    r = mda.analysis.rms.RMSD(u, ref, select=align_range, groupselections=[analysis_range])
    r.run()

    rmsd = r.results.rmsd.T
    rmsd_dict = {'frame': rmsd[1], f'{align_range}': rmsd[2], f"{analysis_range}": rmsd[3]}

    return rmsd_dict


def rmsd_kde(rmsd_data: list, input_dict: dict, widget) -> dict:
    """
    Perform Kernel Density Estimation (KDE) on RMSD data and identify the most distant modes.

    Parameters:
    rmsd_data (list): A list of RMSD values.
    input_dict (dict): A dictionary containing job-related metadata (jobname, max_seq, extra_seq, analysis range, etc.).
    widget (object): A widget object for interactive plotting.

    Returns:
    dict: A dictionary containing information about the identified modes, including their indices, values, and densities.
    """
    jobname = input_dict['jobname']
    max_seq = input_dict['max_seq']
    extra_seq = input_dict['extra_seq']
    rmsd_range = input_dict['analysis_range']
    align_range = input_dict['align_range']
    output_path = input_dict['output_path']

    # Calculate KDE
    kde = gaussian_kde(rmsd_data, bw_method='silverman')
    x_vals = np.linspace(min(rmsd_data), max(rmsd_data), 1000)
    kde_vals = kde(x_vals)

    # Find modes (peaks)
    peaks, _ = find_peaks(kde_vals)
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
        array = np.asarray(rmsd_data).flatten()  # Flatten the array to ensure it's 1-dimensional
        target_value = np.asarray(mode).item()  # Ensure target_value is a scalar

        mode_index = (np.abs(array - target_value)).argmin()
        mode_results = {'mode_index': mode_index,
                        'mode_value': mode,
                        'mode_density': mode_density}

        mode_n += 1

        modes_dict[f'mode_{mode_n}'] = mode_results
    
    if widget:
        plot_item = widget.add_plot(x_vals, kde_vals, title=f'{jobname} {max_seq} {extra_seq}', x_label='RMSD (Å)', y_label='Density', label='KDE')
        widget.add_scatter(plot_item, modes, kde_vals[peaks], label='Modes')
        
    for mode in most_distant_modes:
        lines=pg.InfiniteLine(pos=mode, angle=90, pen=pg.mkPen('b', style=QtCore.Qt.DashLine), label='Most Distant Mode')
        plot_item.addItem(lines)
        
    if output_path:
        # Plot KDE and mark modes
        plt.figure(figsize=(6, 3))
        plt.plot(x_vals, kde_vals, label='KDE')
        plt.scatter(modes, kde_vals[peaks], color='red', zorder=5, label='Modes')

        for mode in most_distant_modes:
            plt.axvline(mode, color='blue', linestyle='--', label='Most Distant Mode')

        plt.title(f'{jobname} {max_seq} {extra_seq}', fontsize=16)
        #plt.suptitle(f'{rmsd_range} after alignment to {align_range}', fontsize=10)
        plt.xlabel('RMSD (Å)', fontsize=14)
        plt.ylabel('Density', fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)  # Major ticks
        plt.tick_params(axis='both', which='minor', labelsize=12)  # Minor ticks (if any)
        plt.grid(False)

        # Show the plot
        plt.tight_layout()

        full_output_path = (f"{output_path}/"
                            f"{jobname}/"
                            f"analysis/"
                            f"rmsd_1d/"
                            f"{jobname}_"
                            f"{max_seq}_"
                            f"{extra_seq}_"
                            f"rmsd_1d_kde.png")

        plt.savefig(full_output_path, dpi=300)
        plt.close()

    return modes_dict


def rmsd_mode_analysis(prediction_dicts, input_dict, ref1d, widget=None):
    """
    Perform 1D RMSD mode analysis for each prediction in the provided dictionary.

    Parameters:
    prediction_dicts (dict): A dictionary containing prediction data with associated MDAnalysis Universes.
    input_dict (dict): A dictionary containing job-related metadata (jobname, analysis range, alignment range, etc.).
    ref1d (str): Path to the reference PDB file for RMSD calculation.
    widget (object): A widget object for interactive plotting.

    Returns:
    dict: The updated prediction_dicts with calculated RMSD data and identified modes.
    """
    jobname = input_dict['jobname']
    rmsd_range = input_dict['analysis_range']
    align_range = input_dict['align_range']
    output_path = input_dict['output_path']

    print('')

    with tqdm(total=len(prediction_dicts), bar_format=TQDM_BAR_FORMAT) as pbar:
        for prediction in prediction_dicts:
            pbar.set_description(f'Running 1D RMSD analysis for {prediction}')
            universe = prediction_dicts[prediction]['mda_universe']
            max_seq = prediction_dicts[prediction]['max_seq']
            extra_seq = prediction_dicts[prediction]['extra_seq']
            input_dict['max_seq'] = max_seq
            input_dict['extra_seq'] = extra_seq

            prediction_dicts[prediction]['rmsd_data'] = calculate_rmsd(universe,
                                                                       ref1d,
                                                                       align_range=align_range,
                                                                       analysis_range=rmsd_range)

            rmsd_df = pd.DataFrame.from_dict(prediction_dicts[prediction]['rmsd_data'], orient='columns')

            full_output_path = (f"{output_path}/"
                                f"{jobname}/"
                                f"analysis/"
                                f"rmsd_1d/"
                                f"{jobname}_"
                                f"{max_seq}_"
                                f"{extra_seq}_"
                                f"rmsd_1d_df.csv")

            rmsd_df.to_csv(full_output_path, index=False)

            all_rmsd_modes = {}
            rmsd_modes = rmsd_kde(prediction_dicts[prediction]['rmsd_data'][rmsd_range], input_dict, widget)

            for mode in rmsd_modes:
                mode_index = rmsd_modes[mode]['mode_index']
                mode_value = int(rmsd_modes[mode]['mode_value'])
                frame_to_save = prediction_dicts[prediction]['rmsd_data']['frame'][mode_index]
                universe.trajectory[int(frame_to_save)]

                full_pdb_path = (f"{output_path}/"
                                 f"{jobname}/"
                                 f"analysis/"
                                 f"representative_structures/"
                                 f"rmsd_1d/"
                                 f"{jobname}_"
                                 f'{max_seq}_'
                                 f'{extra_seq}_'
                                 f'frame_'
                                 f'{mode_index}_'
                                 f'rmsd_vs_ref_'
                                 f'{mode_value}A.pdb')

                with mda.Writer(full_pdb_path, universe.atoms.n_atoms) as W:
                    W.write(universe.atoms)

                rmsd_modes[mode]['pdb_filename'] = full_pdb_path

                all_rmsd_modes[rmsd_range] = rmsd_modes

            prediction_dicts[prediction]['rmsd_data']['rmsd_modes'] = all_rmsd_modes
            pbar.update(n=1)

    return prediction_dicts


def build_dataset_rmsd_modes(results_dict, input_dict):
    """
    Build a dataset from the RMSD mode analysis results and save it as a CSV file.

    Parameters:
    results_dict (dict): A dictionary containing the results of the RMSD mode analysis.
    input_dict (dict): A dictionary containing job-related metadata (jobname, analysis range, output path, etc.).

    Returns:
    None: The function saves the dataset as a CSV file in the specified output directory.
    """
    jobname = input_dict['jobname']
    rmsd_range = input_dict['analysis_range']
    output_path = input_dict['output_path']

    trials = []
    rmsd_range_labels = []
    mode_labels = []
    mode_indexes = []
    mode_densities = []
    mode_values = []
    pdb_filenames = []

    for result in results_dict:
        mode_data = results_dict[result]['rmsd_data']['rmsd_modes'][rmsd_range]
        for mode in mode_data:
            trials.append(result)
            rmsd_range_labels.append(rmsd_range)
            mode_labels.append(mode)
            mode_indexes.append(results_dict[result]['rmsd_data']['rmsd_modes'][rmsd_range][mode]['mode_index'])
            mode_densities.append(results_dict[result]['rmsd_data']['rmsd_modes'][rmsd_range][mode]['mode_density'])
            mode_values.append(results_dict[result]['rmsd_data']['rmsd_modes'][rmsd_range][mode]['mode_value'])
            pdb_filenames.append(results_dict[result]['rmsd_data']['rmsd_modes'][rmsd_range][mode]['pdb_filename'])

    df = pd.DataFrame({'mode_label': mode_labels,
                       'representative_Frame': mode_indexes,
                       'RMSD': mode_values,
                       'mode_peak_density': mode_densities,
                       'trial': trials,
                       'rmsd_range_label': rmsd_range_labels,
                       'pdb_filename': pdb_filenames
                       })

    full_output_path = (f"{output_path}/"
                        f"{jobname}/"
                        f"analysis/"
                        f"rmsd_1d/"
                        f"{jobname}_"
                        f"rmsd_1d_analysis_results.csv")

    print(f"\nSaving {jobname} RMSD_1D analysis results to {full_output_path}\n")

    print(f"Representative Structures for {jobname} RMSD 1D analysis saved at {output_path}/{jobname}/analysis/"
          f"representative_structures/rmsd_1d/\n")

    df.to_csv(full_output_path, index=False)
