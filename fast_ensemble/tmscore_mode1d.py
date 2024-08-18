import os
import warnings
from fast_ensemble.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config
from fast_ensemble.ensemble_analysis.tmscore import build_dataset_tmscore_modes, tmscore_mode_analysis

warnings.filterwarnings("ignore")
import argparse

def run_tmscore_analysis(config, plot_widget=None):
    """
    Run TM-Score analysis based on the provided configuration.

    Parameters:
    config (dict): Configuration dictionary containing parameters for the analysis.
    plot_widget (object, optional): Widget for displaying results (default is None).

    Raises:
    NotADirectoryError: If the specified output path is not a directory.
    """

    # Retrieve configuration values
    output_path = config.get('output_path')
    predictions_path = config.get('predictions_path')
    seq_pairs = config.get('seq_pairs')
    jobname = config.get('jobname')
    slice_predictions = config.get('slice_predictions')
    ref1 = config.get('ref1')
    engine = config.get('engine')
    starting_residue = config.get('starting_residue')

    # Check if the output path is a valid directory
    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    # Set default predictions path if not provided
    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    # Create necessary directories
    create_directory(f'{output_path}/{jobname}/analysis/tmscore_1d')
    create_directory(f'{output_path}/{jobname}/analysis/representative_structures/tmscore_1d')

    # Display configurations
    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Predictions Path: {predictions_path}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")
    print(f"Engine: {engine}")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    if ref1:
        print(f"Reference: {ref1}")
    if slice_predictions:
        print(f"Setting Analysis Range to: {slice_predictions}")
    print("***************************************************************\n")

    # Prepare input dictionary
    input_dict = {
        'jobname': jobname,
        'predictions_path': predictions_path,
        'output_path': output_path,
        'seq_pairs': seq_pairs
    }

    # Load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)
    print("Load predictions to RAM")

    # Run 1D TM-Score analysis
    tmscore_mode_analysis_dict = tmscore_mode_analysis(pre_analysis_dict, input_dict, ref1, slice_predictions, plot_widget)
    print("Run 1D TM-Score analysis")

    # Build results dataset and save to disk
    build_dataset_tmscore_modes(tmscore_mode_analysis_dict, input_dict)


def main():
    """
    Main function to parse arguments and run TM-Score analysis.
    """

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Run TM-Score analysis for a set of predictions.")

    parser.add_argument('--config_file', type=str, default='config.json',
                        help="OPTIONAL: Path to load configuration from file (default: config.json)")
    parser.add_argument('--output_path', type=str, help="Path to save results to")
    parser.add_argument('--predictions_path', type=str, help="Path to read predictions from")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--seq_pairs', type=str, help="A list of [max_seq, extra_seq] pairs used for predictions")
    parser.add_argument('--starting_residue', type=int, help="Sets the starting residue for reindexing (predictions are usually 1-indexed)")
    parser.add_argument('--slice_predictions', type=str, help="The slice range of predictions to analyze")
    parser.add_argument('--ref1', type=str, help="Reference structure for TM-Score calculations")
    parser.add_argument('--engine', type=str, help="Engine used to generate predictions (e.g., AlphaFold2, OpenFold)")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    config.update({k: v for k, v in vars(args).items() if v is not None})

    # Run TM-Score analysis with the provided configuration
    run_tmscore_analysis(config)


if __name__ == "__main__":
    main()
