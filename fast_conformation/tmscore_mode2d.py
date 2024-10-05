import os
import warnings
from fast_conformation.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config, auto_select_2d_references
from fast_conformation.ensemble_analysis.twotmscore import TwoTMScore
import argparse
warnings.filterwarnings("ignore")


def run_2d_tmscore_analysis(config, widget=None):
    """
    Run 2D TM-Score analysis based on the provided configuration.

    Parameters:
    config (dict): Configuration dictionary containing parameters for the analysis.
    widget (object, optional): Widget for displaying results (default is None).

    Raises:
    NotADirectoryError: If the specified output path is not a directory.
    """

    # Retrieve configuration values
    output_path = config.get('output_path')
    predictions_path = config.get('predictions_path')
    mode_results = config.get('mode_results')
    seq_pairs = config.get('seq_pairs')
    jobname = config.get('jobname')
    ref2d1 = config.get('ref2d1')
    ref2d2 = config.get('ref2d2')
    starting_residue = config.get('starting_residue')
    slice_predictions = config.get('slice_predictions')
    engine = config.get('engine')
    n_stdevs = config.get('n_stdevs')
    n_clusters = config.get('n_clusters')

    # Check if the output path is a valid directory
    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    # Set default predictions path if not provided
    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    # Set default mode results path if not provided
    if not mode_results:
        mode_results = f'{output_path}/{jobname}/analysis/tmscore_1d/{jobname}_tmscore_1d_analysis_results.csv'

    # Auto-select references if not provided
    if not ref2d1 and not ref2d2:
        ref2d1, ref2d2 = auto_select_2d_references(mode_results, 'tmscore')

    # Create necessary directories
    create_directory(f'{output_path}/{jobname}/analysis/tmscore_2d')

    # Display configurations
    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Predictions Path: {predictions_path}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")
    print(f"Engine: {engine}")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    if slice_predictions:
        print(f"Setting Analysis Range to: {slice_predictions}")
    print(f"Reference 1: {ref2d1}")
    print(f"Reference 2: {ref2d2}")
    print(f"Number of Standard Devs. to Consider Point Closeness: {n_stdevs}")
    if n_clusters:
        print(f"Number of Clusters: {n_clusters}")
    else:
        print(f"Number of Clusters: Number of Detected 1D TM-Score Modes + 1")
    print("***************************************************************\n")

    # Prepare input dictionary
    input_dict = {
        'jobname': jobname,
        'output_path': output_path,
        'seq_pairs': seq_pairs,
        'predictions_path': predictions_path
    }

    # Load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # Run 2D TM-Score analysis
    twod = TwoTMScore(pre_analysis_dict, input_dict, widget, ref2d1, ref2d2, slice_predictions)

    # Build results dataset and save to disk
    twod.get_2d_tmscore(mode_results, n_stdevs, n_clusters, output_path)


def main():
    """
    Main function to parse arguments and run 2D TM-Score analysis.
    """

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Run 2D TM-Score analysis for a set of predictions.")

    parser.add_argument('--config_file', type=str, default='config.json',
                        help="OPTIONAL: Path to load configuration from file (default: config.json)")
    parser.add_argument('--output_path', type=str, help="Path to save results to")
    parser.add_argument('--predictions_path', type=str, help="Path to read predictions from")
    parser.add_argument('--mode_results', type=str, help="Path to the mode results CSV file")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--seq_pairs', type=str, help="A list of [max_seq, extra_seq] pairs used for predictions")
    parser.add_argument('--starting_residue', type=int, help="Sets the starting residue for reindexing (predictions are usually 1-indexed)")
    parser.add_argument('--slice_predictions', type=str, help="The slice range of predictions to analyze")
    parser.add_argument('--engine', type=str, help="Engine used to generate predictions (e.g., AlphaFold2, OpenFold)")
    parser.add_argument('--ref2d1', type=str, help="First reference structure for TM-Score calculations")
    parser.add_argument('--ref2d2', type=str, help="Second reference structure for TM-Score calculations")
    parser.add_argument('--n_stdevs', type=int, help="Number of standard deviations to consider when calculating close points to fit curve")
    parser.add_argument('--n_clusters', type=int, help="Number of clusters to consider for TM-Score analysis")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    config.update({k: v for k, v in vars(args).items() if v is not None})

    # Run 2D TM-Score analysis with the provided configuration
    run_2d_tmscore_analysis(config)


if __name__ == "__main__":
    main()
