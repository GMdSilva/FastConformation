import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

import MDAnalysis as mda
from MDAnalysis.analysis import rms


def calculate_rmsd(u: mda.Universe,
                   ref: mda.Universe = None,
                   align_range: str = "backbone",
                   analysis_range: str = "backbone") -> dict:
    """
    Calculate the rmsd of a universe
    Parameters:
        u (mda.Universe): Universe to be analyzed
        ref (mda.Universe): Reference universe to be analyzed
        align_range (str):
        analysis_range (str):
    Returns:
        rmsd_dict (dict):
    """

    # if reference is not supplied, use first frame
    if ref is None:
        ref = u.select_atoms('protein')

    print(f"Calculating RMSD of '{analysis_range}' after '{align_range}' alignment.")

    R = mda.analysis.rms.RMSD(u, ref, select=align_range, groupselections=[analysis_range])
    R.run()

    rmsd = R.results.rmsd.T
    rmsd_dict = {'frame': rmsd[1], f'{align_range}': rmsd[2]}
    rmsd_dict[analysis_range] = rmsd[3]

    return rmsd_dict


def rmsd_kde(rmsd_data: list, input_dict: dict) -> dict:
    jobname = input_dict['jobname']
    max_seq = input_dict['max_seq']
    extra_seq = input_dict['extra_seq']
    rmsd_range = input_dict['rmsd_range']
    align_range = input_dict['align_range']
    output_path = input_dict['output_path']
    rmsd_range_n = input_dict['rmsd_range_name']

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

    print(f"Plotting RMSD mode data.")

    # Plot KDE and mark modes
    plt.figure(figsize=(6, 3))
    plt.plot(x_vals, kde_vals, label='KDE')
    plt.scatter(modes, kde_vals[peaks], color='red', zorder=5, label='Modes')

    for mode in most_distant_modes:
        plt.axvline(mode, color='blue', linestyle='--', label='Most Distant Mode')

    plt.title(f'{jobname} max_seq: {max_seq} extra_seq: {extra_seq}', fontsize=16)
    plt.suptitle(f'{rmsd_range} after alignment by {align_range}')
    plt.xlabel('RMSD ($\AA$)', fontsize=14)
    plt.ylabel('Density', fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)  # Major ticks
    plt.tick_params(axis='both', which='minor', labelsize=12)  # Minor ticks (if any)
    plt.grid(False)

    # Show the plot
    plt.tight_layout()
    plt.savefig(f"{output_path}/plots/{jobname}_{max_seq}_{extra_seq}_rmsd_kde_{rmsd_range_n}.png", dpi=300)
    print(f"Saved RMSD mode plot to {output_path}/plots/{jobname}_{max_seq}_{extra_seq}_rmsd_kde_{rmsd_range_n}.png")
    plt.close()

    return modes_dict


def rmsd_mode_analysis(prediction_dicts, input_dict):
    jobname = input_dict['jobname']
    rmsd_range = input_dict['rmsd_range']
    align_range = input_dict['align_range']
    output_path = input_dict['output_path']

    for prediction in prediction_dicts:
        print(f"Finding representative structures for {prediction}.")

        universe = prediction_dicts[prediction]['mda_universe']
        max_seq = prediction_dicts[prediction]['max_seq']
        extra_seq = prediction_dicts[prediction]['extra_seq']
        input_dict['max_seq'] = max_seq
        input_dict['extra_seq'] = extra_seq

        prediction_dicts[prediction]['rmsd_data'] = calculate_rmsd(universe,
                                                                   align_range=align_range,
                                                                   analysis_range=rmsd_range)

        all_rmsd_modes = {}
        rmsd_modes = rmsd_kde(prediction_dicts[prediction]['rmsd_data'][rmsd_range], input_dict)

        for mode in rmsd_modes:
            mode_index = rmsd_modes[mode]['mode_index']
            mode_value = int(rmsd_modes[mode]['mode_value'])
            frame_to_save = prediction_dicts[prediction]['rmsd_data']['frame'][mode_index]
            universe.trajectory[int(frame_to_save)]
            # Save the frame to a file
            pdb_filename = (f'{output_path}/representative_structures/'
                            f'{jobname}_{max_seq}_{extra_seq}_frame_{mode_index}_rmsd_vs_ref_{mode_value}A.pdb')

            print(f"Saving {mode} representative structure to "
                  f"{output_path}/"
                  f"representative_structures/"
                  f"{jobname}_{max_seq}_{extra_seq}_frame_{mode_index}_rmsd_vs_ref_{mode_value}A.pdb")

            with mda.Writer(pdb_filename, universe.atoms.n_atoms) as W:
                W.write(universe.atoms)

            rmsd_modes[mode]['pdb_filename'] = pdb_filename

            all_rmsd_modes[rmsd_range] = rmsd_modes

        prediction_dicts[prediction]['rmsd_data']['rmsd_modes'] = all_rmsd_modes

    return prediction_dicts


def build_dataset_rmsd_modes(results_dict, input_dict):
    jobname = input_dict['jobname']
    rmsd_range = input_dict['rmsd_range']
    align_range = input_dict['align_range']
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
                       'RMSD_value_(A)': mode_values,
                       'mode_density': mode_densities,
                       'trial': trials,
                       'rmsd_range_label': rmsd_range_labels,
                       'pdb_filename': pdb_filenames
                       })

    print(f"Saving {jobname} CSV data structure to {output_path}/datasets/{jobname}_rmsd_mode_analysis_results.csv")

    df.to_csv(f"{output_path}/datasets/{jobname}_rmsd_mode_analysis_results.csv")

