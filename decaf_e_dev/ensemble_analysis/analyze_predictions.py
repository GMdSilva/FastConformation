import os
import json
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
import shutil
import warnings
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

from sklearn.cluster import KMeans
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

warnings.filterwarnings('ignore')


def create_directory(path):
    # Check if the directory already exists
    if os.path.exists(path):
        # If it exists, delete it
        shutil.rmtree(path)
    # Create the new directory
    os.makedirs(path)
    print(f"Directory '{path}' created successfully.")


def load_config(config_file):
    try:
        with open(config_file, 'r') as file:
            config = json.load(file)
        return config
    except FileNotFoundError:
        print(f"Configuration file {config_file} not found.")
    except json.JSONDecodeError:
        print(f"Error decoding JSON from the configuration file {config_file}.")


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
    print(f"Loaded Universe with {len(u.trajectory)} frames")
    return u


def load_predictions(output_path, seq_pairs, jobname):
    pred_dict = {}
    for seq_pair in seq_pairs:
        max_seq = seq_pair[0]
        extra_seq = seq_pair[1]
        folder_path = f"{output_path}/{jobname}_{max_seq}_{extra_seq}"

        print(f"Loading predictions at {folder_path}")
        universe = load_pdb_files_as_universe(folder_path)

        params = {'max_seq': max_seq,
                  'extra_seq': extra_seq,
                  'jobname': jobname,
                  'mda_universe': universe}
        pred_dict[f'{jobname}_{max_seq}_{extra_seq}'] = params
    return pred_dict


def calculate_rmsd(u, ref=None, align_range="backbone", rmsd_ranges=("backbone and resid 150-170")):
    if ref is None:
        ref = u.select_atoms('protein')

    R = mda.analysis.rms.RMSD(u, ref, select=align_range, groupselections=rmsd_ranges)
    R.run()

    rmsd = R.results.rmsd.T

    rmsd_dict = {'frame': rmsd[1],
                 f'{align_range}': rmsd[2]}

    for rmsd_values, rmsd_range in zip(rmsd[3:], rmsd_ranges):
        rmsd_dict[rmsd_range] = rmsd_values

    return rmsd_dict


def rmsd_kde(rmsd_data, jobname, max_seq, extra_seq, rmsd_range, output_path, rmsd_range_n):
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

    # Find the two most populated modes
    sorted_modes = sorted(zip(modes, kde_vals[peaks]), key=lambda x: x[1], reverse=True)
    most_populated_modes = np.array([mode[0] for mode in sorted_modes[:2]])

    # Plot KDE and mark modes
    plt.figure(figsize=(8, 3))
    plt.plot(x_vals, kde_vals, label='KDE')
    plt.scatter(modes, kde_vals[peaks], color='red', zorder=5, label='Modes')

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

    for mode in most_distant_modes:
        plt.axvline(mode, color='blue', linestyle='--', label='Most Distant Mode')

    # Customize the plot
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
    plt.close()

    return modes_dict


def rmsd_mode_analysis(prediction_dicts, align_range, rmsd_ranges, output_path):
    for prediction in prediction_dicts:
        universe = prediction_dicts[prediction]['mda_universe']
        prediction_dicts[prediction]['rmsd_data'] = calculate_rmsd(universe,
                                                                   align_range=align_range,
                                                                   rmsd_ranges=rmsd_ranges)

        all_rmsd_modes = {}
        for rmsd_range_n, rmsd_range in enumerate(rmsd_ranges):
            jobname = prediction_dicts[prediction]['jobname']
            max_seq = prediction_dicts[prediction]['max_seq']
            extra_seq = prediction_dicts[prediction]['extra_seq']

            rmsd_modes = rmsd_kde(prediction_dicts[prediction]['rmsd_data'][rmsd_range],
                                          jobname,
                                          max_seq,
                                          extra_seq,
                                          rmsd_range,
                                          output_path,
                                          rmsd_range_n)


            for mode in rmsd_modes:
                mode_index = rmsd_modes[mode]['mode_index']
                mode_value = int(rmsd_modes[mode]['mode_value'])
                frame_to_save = prediction_dicts[prediction]['rmsd_data']['frame'][mode_index]
                universe.trajectory[int(frame_to_save)]
                # Save the frame to a file
                pdb_filename = (f'{output_path}/representative_structures/'
                                f'{jobname}_{max_seq}_{extra_seq}_frame_{mode_index}_rmsd_vs_ref_{mode_value}A.pdb')

                with mda.Writer(pdb_filename, universe.atoms.n_atoms) as W:
                    W.write(universe.atoms)
                rmsd_modes[mode]['pdb_filename'] = pdb_filename

            all_rmsd_modes[rmsd_range] = rmsd_modes

        prediction_dicts[prediction]['rmsd_data']['rmsd_modes'] = all_rmsd_modes

    return prediction_dicts


def calculate_rmsf_and_call_peaks(jobname,
                   prediction_dicts,
                   align_range,
                   output_path):

    for result in prediction_dicts:
        u = prediction_dicts[result]['mda_universe']
        max_seq = prediction_dicts[result]['max_seq']
        extra_seq = prediction_dicts[result]['extra_seq']

        average = align.AverageStructure(u, u, select=align_range, ref_frame=0).run()
        ref = average.results.universe
        aligner = align.AlignTraj(u, ref, select=align_range, in_memory=True).run()
        atom_sel = u.select_atoms(align_range)
        R = rms.RMSF(atom_sel).run()

        rmsf_values = R.results.rmsf
        resids = atom_sel.resids

        mean_rmsf = np.mean(rmsf_values)
        std_rmsf = np.std(rmsf_values)
        threshold = mean_rmsf + 2 * std_rmsf  # You can adjust the multiplier as needed

        # Detect peaks
        peaks, properties = find_peaks(rmsf_values, width=3, prominence=1, height=threshold)

        plt.figure(figsize=(8, 3))
        plt.title(f'{jobname} {max_seq} {extra_seq} aligned to {align_range}', fontsize=16)
        plt.xlabel('Data Point (Not Residue Number)', fontsize=14)
        plt.ylabel('RMSF ($\AA$)', fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)  # Major ticks
        plt.tick_params(axis='both', which='minor', labelsize=12)  # Minor ticks (if any)

        plt.plot(rmsf_values)
        plt.plot(peaks, rmsf_values[peaks], "x")
        plt.vlines(x=peaks, ymin=rmsf_values[peaks] - properties["prominences"],
                   ymax=rmsf_values[peaks], color="C1")
        plt.hlines(y=properties["width_heights"], xmin=properties["left_ips"],
                   xmax=properties["right_ips"], color="C1")

        plt.tight_layout()
        plt.savefig(f'{output_path}/plots/{jobname}_{max_seq}_{extra_seq}_rmsf_{align_range}_peaks.png', dpi=300)

        detected_peaks = {}
        peak_counter = 0
        # Extract resid values for each peak
        for i, peak in enumerate(peaks):
            left_ip = int(properties["left_ips"][i])
            right_ip = int(properties["right_ips"][i])
            peak_resids = resids[left_ip:right_ip + 1]

            peak_prop = {f'starting_residue': peak_resids[0],
                         f'ending_residue': peak_resids[-1],
                         f'length': len(peak_resids),
                         f'peak_value': rmsf_values[peak],
                         f'prominence': properties["prominences"][i],
                         f'width_height': properties["width_heights"][i]}

            peak_counter += 1

            detected_peaks[f'detected_peak_{peak_counter}'] = peak_prop

        prediction_dicts[result]['detected_peaks'] = detected_peaks

    return prediction_dicts


def calculate_rmsf_multiple(jobname,
                   prediction_dicts,
                   align_range,
                   output_path,
                   colors=['red', 'blue', 'green', 'purple', 'orange']):
    labels = []
    plt.figure(figsize=(8, 3))
    plt.title(f'{jobname} aligned to {align_range}', fontsize=16)
    plt.xlabel('Residue number', fontsize=14)
    plt.ylabel('RMSF ($\AA$)', fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)  # Major ticks
    plt.tick_params(axis='both', which='minor', labelsize=12)  # Minor ticks (if any)

    for idx, (result, data) in enumerate(prediction_dicts.items()):
        u = data['mda_universe']
        average = align.AverageStructure(u, u, select=align_range, ref_frame=0).run()
        ref = average.results.universe

        aligner = align.AlignTraj(u, ref, select=align_range, in_memory=True).run()

        atom_sel = u.select_atoms(align_range)
        R = rms.RMSF(atom_sel).run()

        plt.plot(atom_sel.resids, R.results.rmsf, color=colors[idx % len(colors)], label=result)
        labels.append(result)  # Add the result to labels list

    plt.legend()  # Add the legend at the end
    plt.tight_layout()
    plt.savefig(f"{output_path}/plots/{jobname}_rmsf_{align_range}.png", dpi=300)
    plt.close()


def build_dataset_rmsd_modes(results_dict, rmsd_ranges, output_path):
    trials = []
    rmsd_range_labels = []
    mode_labels = []
    mode_indexes = []
    mode_densities = []
    mode_values = []
    pdb_filenames = []

    for result in results_dict:
        for rmsd_range in rmsd_ranges:
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

    df.to_csv(f"{output_path}/datasets/{jobname}_rmsd_mode_analysis_results.csv")


def build_dataset_rmsf_peaks(results_dict, output_path):
    trials = []
    peak_labels = []
    peak_starting_residues = []
    peak_ending_residues = []
    peak_lengths = []
    peak_rmsf_values = []
    peak_prominences = []
    peak_width_heights = []

    for result in results_dict:
        peak_data = results_dict[result]['detected_peaks']
        for peak in peak_data:
            trials.append(result)
            peak_labels.append(peak)
            peak_starting_residues.append(results_dict[result]['detected_peaks'][peak]['starting_residue'])
            peak_ending_residues.append(results_dict[result]['detected_peaks'][peak]['ending_residue'])
            peak_rmsf_values.append(results_dict[result]['detected_peaks'][peak]['peak_value'])
            peak_prominences.append(results_dict[result]['detected_peaks'][peak]['prominence'])
            peak_lengths.append(results_dict[result]['detected_peaks'][peak]['length'])
            peak_width_heights.append(results_dict[result]['detected_peaks'][peak]['width_height'])

    df = pd.DataFrame({'peak_label': peak_labels,
                       'starting_residue': peak_starting_residues,
                       'ending_residue': peak_ending_residues,
                       'length': peak_lengths,
                       'rmsf_value': peak_rmsf_values,
                       'prominence': peak_prominences,
                       'width_height': peak_width_heights,
                       'trial': trials})

    df.to_csv(f"{output_path}/datasets/{jobname}_rmsf_peak_calling_results.csv")


def half_parabola(x, a, b, c):
    return a * x**2 + b * x + c


def get_2d_rmsd(rmsd_mode_df_path, prediction_dicts):
    df = pd.read_csv(rmsd_mode_df_path)

    unique_trials = df['trial'].unique()

    # Iterate and subset by trial
    for trial in unique_trials:
        trial_subset = df[df['trial'] == trial]

        rmsd_range = trial_subset['rmsd_range_label'].iloc[-1]
        n_modes = len(trial_subset['mode_label'])
        ref_gr = trial_subset['pdb_filename'].iloc[0]
        ref_alt = trial_subset['pdb_filename'].iloc[-1]

        universe = prediction_dicts[trial]['mda_universe']
        max_seq = prediction_dicts[trial]['max_seq']
        extra_seq = prediction_dicts[trial]['extra_seq']
        jobname = prediction_dicts[trial]['jobname']

        ref_gr_u = mda.Universe(ref_gr)
        ref_alt_u = mda.Universe(ref_alt)

        rmsd_gr = calculate_rmsd(universe, ref=ref_gr_u, align_range='backbone', rmsd_ranges=[rmsd_range])
        rmsd_alt = calculate_rmsd(universe, ref=ref_alt_u, align_range='backbone', rmsd_ranges=[rmsd_range])

        rmsd_2d_data = np.array([rmsd_gr[rmsd_range], rmsd_alt[rmsd_range]])
        rmsd_2d_data = rmsd_2d_data.T

        kmeans = KMeans(n_clusters=2)
        kmeans.fit(rmsd_2d_data)

        labels = kmeans.labels_
        centroids = kmeans.cluster_centers_

        if centroids[0, 0] < centroids[1, 0]:
            # If centroid 0 is closer to the origin on the x-axis
            correct_labels = labels
        else:
            # Swap the labels if centroid 1 is closer to the origin on the x-axis
            correct_labels = 1 - labels

        plt.figure(figsize=(5, 4))
        colors = ['r', 'g']  # Extend this if you have more than 6 clusters

        for i in range(2):
            cluster_points = rmsd_2d_data[correct_labels == i]
            plt.scatter(cluster_points[:, 0], cluster_points[:, 1], c=colors[i], label=f'Cluster {i + 1}', alpha=0.6)

        if len(rmsd_2d_data) > 0:
            popt, _ = curve_fit(half_parabola, rmsd_2d_data[:, 0], rmsd_2d_data[:, 1])
            fit_x = np.linspace(min(rmsd_2d_data[:, 0]), max(rmsd_2d_data[:, 0]), 100)
            fit_y = half_parabola(fit_x, *popt)
            plt.plot(fit_x, fit_y, c='b', linestyle='--', label='Overall Fit')

            # Calculate goodness of fit metrics
            y_true = rmsd_2d_data[:, 1]
            y_pred = half_parabola(rmsd_2d_data[:, 0], *popt)
            r2 = r2_score(y_true, y_pred)
            rmse = np.sqrt(mean_squared_error(y_true, y_pred))
            mae = mean_absolute_error(y_true, y_pred)

            print(f'Goodness of fit metrics for trial {trial}:')
            print(f'R²: {r2:.4f}')
            print(f'RMSE: {rmse:.4f}')
            print(f'MAE: {mae:.4f}')

        plt.scatter(centroids[:, 0], centroids[:, 1], s=100, c='k', marker='X', label='Centroids')

        # Adding title and labels
        plt.title(f'{jobname} {max_seq} {extra_seq} "{rmsd_range}"', fontsize=16)
        plt.xlabel('RMSD vs. Ref1 (A)', fontsize=14)
        plt.ylabel('RMSD vs. Ref2 (A)', fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)  # Major ticks
        plt.tick_params(axis='both', which='minor', labelsize=12)  # Minor ticks (if any)
        plt.tight_layout()
        plt.show()

        unique_labels = set(correct_labels)

        cluster_counts = {i: 0 for i in unique_labels}

        for k, col in zip(unique_labels, colors):
            class_member_mask = (correct_labels == k)
            xy = rmsd_2d_data[class_member_mask]
            cluster_counts[k] = len(xy)

        # Print cluster populations
        for cluster_id, count in cluster_counts.items():
            print(f'Cluster {cluster_id}: {count} points')



config = load_config('config.json')

# Override config with command line arguments if provided
output_path = config.get('output_path')
seq_pairs = config.get('seq_pairs')
jobname = config.get('jobname')
align_range = 'backbone'
config_file = "config.json"
#
# create_directory(f'{output_path}/representative_structures')
# create_directory(f'{output_path}/plots')
# create_directory(f'{output_path}/datasets')

pre_analysis_dict = load_predictions(output_path, seq_pairs, jobname)
#
# calculate_rmsf_multiple(jobname, pre_analysis_dict, align_range, output_path)
# rmsf_peak_calling_dict = calculate_rmsf_and_call_peaks(jobname, pre_analysis_dict, align_range, output_path)
# build_dataset_rmsf_peaks(rmsf_peak_calling_dict, output_path)
#
#
# for result in rmsf_peak_calling_dict:
#     for peak in rmsf_peak_calling_dict[result]['detected_peaks']:
#         starting_residue = rmsf_peak_calling_dict[result]['detected_peaks'][peak]['starting_residue']
#         ending_residue = rmsf_peak_calling_dict[result]['detected_peaks'][peak]['ending_residue']
#
#
# rmsd_ranges = [f'backbone and resid {starting_residue}-{ending_residue}']
#
#
# rmsd_mode_analysis_dict = rmsd_mode_analysis(pre_analysis_dict, align_range, rmsd_ranges, output_path)
# build_dataset_rmsd_modes(rmsd_mode_analysis_dict, rmsd_ranges, output_path)

get_2d_rmsd('sample_predictions/datasets/abl_wt_rmsd_mode_analysis_results.csv', pre_analysis_dict)
