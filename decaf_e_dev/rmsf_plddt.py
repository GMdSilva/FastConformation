import os
import warnings
from decaf_e_dev.ensemble_analysis.analysis_utils import (
    load_config, create_directory, load_predictions, load_predictions_json
)
from decaf_e_dev.ensemble_analysis.rmsf import (
    calculate_rmsf_multiple, calculate_rmsf_and_call_peaks, build_dataset_rmsf_peaks, plot_plddt_line, plot_plddt_rmsf_corr
)

warnings.filterwarnings("ignore")

def run_rmsf_analysis(config):

    jobname = config.get('jobname')
    output_path = config.get('output_path')
    seq_pairs = config.get('seq_pairs')
    predictions_path = config.get('predictions_path')
    engine = config.get('engine')
    align_range = config.get('align_range')
    detect_mobile = config.get('detect_mobile')
    peak_width = config.get('peak_width')
    peak_prominence = config.get('peak_prominence')
    peak_height = config.get('peak_height')
    starting_residue = config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    create_directory(f'{output_path}/{jobname}/analysis/rmsf_plddt')

    if detect_mobile:
        create_directory(f'{output_path}/{jobname}/analysis/mobile_detection')

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Job Name: {jobname}")
    print(f"Output Path: {output_path}")
    print(f"max_seq:extra_seq Pairs: {seq_pairs}")
    print(f"Predictions Path: {predictions_path}")
    print(f"Engine: {engine}")
    print(f"Align Range: {align_range}")
    print(f"Detect Mobile Segments? {detect_mobile}")
    if detect_mobile:
        print(f"Peak Width: {peak_width}")
        print(f"Peak Prominence: {peak_prominence}")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    print("***************************************************************\n")

    # load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # runs RMSF analysis for all predictions
    calculate_rmsf_multiple(jobname, pre_analysis_dict, align_range, output_path)

    # gets pLDDTs for all predictions
    plddt_dict = load_predictions_json(predictions_path, seq_pairs, jobname)
    plot_plddt_line(jobname, plddt_dict, output_path, starting_residue)

    # runs pLDDT/RMSF correlation for all predictions
    plot_plddt_rmsf_corr(jobname, pre_analysis_dict, plddt_dict, output_path)

    # runs mobile residue range detection with RMSF data and save dataset to disk
    if detect_mobile:
        print(f'\nAutomatically detecting mobile segments')
        rmsf_peak_calling_dict = calculate_rmsf_and_call_peaks(jobname, pre_analysis_dict, align_range, output_path, peak_width, peak_prominence, peak_height)

        build_dataset_rmsf_peaks(jobname, rmsf_peak_calling_dict, output_path, engine)
