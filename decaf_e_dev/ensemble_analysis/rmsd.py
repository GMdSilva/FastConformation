import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

import MDAnalysis as mda
from MDAnalysis.analysis import rms

from tqdm import tqdm


TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'


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
    else:
        ref_u = mda.Universe(ref)
        ref = ref_u.select_atoms('protein')

    r = mda.analysis.rms.RMSD(u, ref, select=align_range, groupselections=[analysis_range])
    r.run()

    rmsd = r.results.rmsd.T
    rmsd_dict = {'frame': rmsd[1], f'{align_range}': rmsd[2], f"{analysis_range}": rmsd[3]}

    return rmsd_dict


# Function to create RMSD KDE and plot
def rmsd_kde(rmsd_data, input_dict, plot_widget):
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

    # Plot KDE and mark modes
    plot_widget.plot(x_vals, kde_vals, f'{jobname} max_seq: {max_seq} extra_seq: {extra_seq}', 'RMSD ($\AA$)', 'Density')

    # Save plot
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

    # Return modes (for further analysis if needed)
    return modes

def rmsd_mode_analysis(prediction_dicts, input_dict, ref1d):
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
            rmsd_modes = rmsd_kde(prediction_dicts[prediction]['rmsd_data'][rmsd_range], input_dict)

            for mode in rmsd_modes:
                mode_index = rmsd_modes[mode]['mode_index']
                mode_value = int(rmsd_modes[mode]['mode_value'])
                frame_to_save = prediction_dicts[prediction]['rmsd_data']['frame'][mode_index]
                universe.trajectory[int(frame_to_save)]
                # Save the frame to a file

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
