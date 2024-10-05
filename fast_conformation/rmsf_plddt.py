import os
import warnings
from fast_conformation.ensemble_analysis.analysis_utils import (
    load_config, create_directory, load_predictions, load_predictions_json
)
from fast_conformation.ensemble_analysis.rmsf import (
    calculate_rmsf_multiple, calculate_rmsf_and_call_peaks, build_dataset_rmsf_peaks, plot_plddt_line, plot_plddt_rmsf_corr
)
import argparse
warnings.filterwarnings("ignore")


def run_rmsf_analysis(config, widget=None):
    """
    Run RMSF analysis based on the provided configuration.

    Parameters:
    config (dict): Configuration dictionary containing parameters for the analysis.
    widget (object, optional): Widget for displaying results (default is None).

    Raises:
    NotADirectoryError: If the specified output path is not a directory.
    """

    # Retrieve configuration values
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

    # Check if the output path is a valid directory
    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    # Create necessary directories
    create_directory(f'{output_path}/{jobname}/analysis/rmsf_plddt')
    if detect_mobile:
        create_directory(f'{output_path}/{jobname}/analysis/mobile_detection')

    # Set default predictions path if not provided
    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    # Display configurations
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

    # Load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # Run RMSF analysis for all predictions
    calculate_rmsf_multiple(jobname, pre_analysis_dict, align_range, output_path, widget)

    # Load pLDDT data
    plddt_dict = load_predictions_json(predictions_path, seq_pairs, jobname)

    # Plot pLDDT line
    plot_plddt_line(jobname, plddt_dict, output_path, starting_residue, widget)

    # Plot pLDDT/RMSF correlation
    plot_plddt_rmsf_corr(jobname, pre_analysis_dict, plddt_dict, output_path, widget)

    # Detect mobile segments if specified
    if detect_mobile:
        print(f'\nAutomatically detecting mobile segments')
        rmsf_peak_calling_dict = calculate_rmsf_and_call_peaks(jobname, pre_analysis_dict, align_range, output_path, peak_width, peak_prominence, peak_height, widget)

        # Build and save the dataset
        build_dataset_rmsf_peaks(jobname, rmsf_peak_calling_dict, output_path, engine)


def main():
    """
    Main function to parse arguments and run RMSF analysis.
    """

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Runs RMSF/pLDDT Analysis for a set of predictions")

    parser.add_argument('--config_file', type=str, help="Path to the configuration file")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--output_path', type=str, help="Path to save results to")
    parser.add_argument('--seq_pairs', type=str, help="A list of [max_seq, extra_seq] pairs previously used to generate the predictions")
    parser.add_argument('--predictions_path', type=str, help="Path to read PDB files of predictions, expects format: /jobname_maxseq_extraseq/ (if not provided, will search automatically based on other parameters)")
    parser.add_argument('--engine', type=str, choices=['alphafold2', 'openfold'], help="The engine previously used to generate predictions (AlphaFold2 or OpenFold), used to find predictions if predictions_path is not supplied")
    parser.add_argument('--starting_residue', type=int, help="Sets the starting residue for reindexing (predictions are usually 1-indexed)")
    parser.add_argument('--align_range', type=str, help="The atom alignment range for RMSF calculations (MDAnalysis Syntax)")
    parser.add_argument('--detect_mobile', type=bool, help="Pass True to detect mobile residue ranges")
    parser.add_argument('--peak_width', type=int, help="Sets the RMSF peak width threshold for mobile residue range detection")
    parser.add_argument('--peak_prominence', type=int, help="Sets the RMSF peak prominence threshold for mobile residue range detection")
    parser.add_argument('--peak_height', type=int, help="Sets the RMSF peak height threshold for mobile residue range detection")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    config.update({k: v for k, v in vars(args).items() if v is not None})

    # Run RMSF analysis with the provided configuration
    run_rmsf_analysis(config)


if __name__ == "__main__":
    main()
