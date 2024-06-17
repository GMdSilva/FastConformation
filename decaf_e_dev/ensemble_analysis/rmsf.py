import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

from scipy.signal import find_peaks

from MDAnalysis.analysis import rms, align


def calculate_rmsf_and_call_peaks(jobname,
                                  prediction_dicts,
                                  align_range,
                                  output_path):

    for result in prediction_dicts:

        u = prediction_dicts[result]['mda_universe']
        max_seq = prediction_dicts[result]['max_seq']
        extra_seq = prediction_dicts[result]['extra_seq']
        print(f"Running peak detection for RMSF of {jobname} {max_seq} {extra_seq} after alignment to {align_range}")

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
        print(f"Saved peak calling with RMSF plot for {output_path}/{jobname} {max_seq} {extra_seq}"
              f" after alignment to {align_range}"
              f" at {output_path}/plots/{jobname}_{max_seq}_{extra_seq}_rmsf_{align_range}_peaks.png")

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

            print(f"Peaks detected for {jobname} {max_seq} {extra_seq}: {peak_prop}")

        prediction_dicts[result]['detected_peaks'] = detected_peaks

    return prediction_dicts


def calculate_rmsf_multiple(jobname,
                   prediction_dicts,
                   align_range,
                   output_path,
                   colors=('red', 'blue', 'green', 'purple', 'orange', 'grey', 'brown', 'cyan', 'magenta')):

    print(f"Calculating RMSF for all results in {output_path}/{jobname}* after alignment to {align_range}")
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
    plt.savefig(f"{output_path}/plots/{jobname}_rmsf_{align_range}_all_results.png", dpi=300)
    print(f"Saved all RMSFs plot to {output_path}/plots/{jobname}_rmsf_{align_range}.png")
    plt.close()


def plot_plddt_rmsf_corr(jobname,
                         prediction_dicts,
                         plddt_dicts,
                         output_path):

    for result, plddt_dict in zip(prediction_dicts, plddt_dicts):
        max_seq = prediction_dicts[result]['max_seq']
        extra_seq = prediction_dicts[result]['extra_seq']

        plddt_data = plddt_dicts[result]['all_plddts']
        arrays = np.array(list(plddt_data.values()))
        plddt_avg = np.mean(arrays, axis=0)

        u = prediction_dicts[result]['mda_universe']
        average = align.AverageStructure(u, u, select='name CA', ref_frame=0).run()
        ref = average.results.universe
        aligner = align.AlignTraj(u, ref, select='name CA', in_memory=True).run()

        atom_sel = u.select_atoms('name CA')
        R = rms.RMSF(atom_sel).run()

        norm = Normalize(vmin=0, vmax=len(R.results.rmsf))
        cmap = cm.get_cmap('viridis')
        colors = cmap(norm(np.arange(len(R.results.rmsf))))

        plt.scatter(R.results.rmsf, plddt_avg, c=colors)

        plt.title(f'{jobname} {max_seq} {extra_seq}', fontsize=16)
        plt.xlabel('C-Alpha RMSF (A)', fontsize=14)
        plt.ylabel('Average pLDDT', fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)  # Major ticks
        plt.tick_params(axis='both', which='minor', labelsize=12)  # Minor ticks (if any)
        plt.tight_layout()
        plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=plt.gca(), label='Residue #')
        plt.savefig(f"{output_path}/plots/{jobname}_{max_seq}_{extra_seq}_plddt_rmsf_corr.png", dpi=300)
        plt.close()



def plot_plddt_line(jobname,
                    plddt_dict,
                    output_path,
                    colors=('red', 'blue', 'green', 'purple', 'orange', 'grey', 'brown', 'cyan', 'magenta')):
    labels = []
    plt.figure(figsize=(8, 3))
    plt.title(f'{jobname}', fontsize=16)
    plt.xlabel('Residue number', fontsize=14)
    plt.ylabel('pLDDT', fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)  # Major ticks
    plt.tick_params(axis='both', which='minor', labelsize=12)  # Minor ticks (if any)

    for idx, (result, data) in enumerate(plddt_dict.items()):
        plddt_data = data['all_plddts']
        arrays = np.array(list(plddt_data.values()))
        plddt_avg = np.mean(arrays, axis=0)
        plt.plot(plddt_avg, color=colors[idx % len(colors)], label=result)
        labels.append(result)  # Add the result to labels list

    plt.legend()  # Add the legend at the end
    plt.tight_layout()
    plt.savefig(f"{output_path}/plots/{jobname}_plddt.png", dpi=300)
    plt.close()


def build_dataset_rmsf_peaks(jobname, results_dict, output_path):
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
    print(f"Saved peak calling output of all results at {output_path}/{jobname}* "
          f"to {output_path}/datasets/{jobname}_rmsf_peak_calling_results.csv")