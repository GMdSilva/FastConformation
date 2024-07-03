import os
import warnings
from decaf_e_dev.ensemble_analysis.analysis_utils import (
    load_config, create_directory, load_predictions, load_predictions_json
)
from decaf_e_dev.ensemble_analysis.rmsf import (
    calculate_rmsf_multiple, calculate_rmsf_and_call_peaks, build_dataset_rmsf_peaks, plot_plddt_line, plot_plddt_rmsf_corr
)

warnings.filterwarnings("ignore")

def run_rmsf_analysis(config_file=None, jobname=None, output_path=None, seq_pairs=None,
                      predictions_path=None, engine=None, align_range=None,
                      detect_mobile=None, peak_width=None, peak_prominence=None,
                      peak_height=None, starting_residue=None):

    # Load configuration from file if provided
    config_file = config_file if config_file else 'config.json'
    config = load_config(config_file)

    # Override config with function arguments if provided
    jobname = jobname if jobname else config.get('jobname')
    output_path = output_path if output_path else config.get('output_path')
    seq_pairs = seq_pairs if seq_pairs else config.get('seq_pairs')
    predictions_path = predictions_path if predictions_path else config.get('predictions_path')
    engine = engine if engine else config.get('engine')
    align_range = align_range if align_range else config.get('align_range')
    detect_mobile = detect_mobile if detect_mobile else config.get('detect_mobile')
    peak_width = peak_width if peak_width else config.get('peak_width')
    peak_prominence = peak_prominence if peak_prominence else config.get('peak_prominence')
    peak_height = peak_height if peak_height else config.get('peak_height')
    starting_residue = starting_residue if starting_residue else config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    create_directory(f'{output_path}/{jobname}/analysis/rmsf_plddt')

    if detect_mobile:
        create_directory(f'{output_path}/{jobname}/analysis/mobile_detection')

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Used Config File? {config_file is not None}")
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
