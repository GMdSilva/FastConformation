import os
import warnings
import argparse
from fast_ensemble.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config
from fast_ensemble.ensemble_analysis.pca import pca_from_ensemble

warnings.filterwarnings("ignore")


def run_pca_analysis(config, widget=None):
    """
    Run PCA analysis based on the provided configuration.

    Parameters:
    config (dict): Configuration dictionary containing parameters for the analysis.
    widget (object, optional): Widget for displaying results (default is None).

    Raises:
    NotADirectoryError: If the specified output path is not a directory.
    """

    # Retrieve configuration values
    predictions_path = config.get('predictions_path')
    output_path = config.get('output_path')
    seq_pairs = config.get('seq_pairs')
    jobname = config.get('jobname')
    align_range = config.get('align_range')
    analysis_range = config.get('analysis_range')
    analysis_range_name = config.get('analysis_range_name')
    engine = config.get('engine')
    n_pca_clusters = config.get('n_pca_clusters')
    starting_residue = config.get('starting_residue')

    # Check if the output path is a valid directory
    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    # Set default predictions path if not provided
    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    # Create necessary directories
    create_directory(f'{output_path}/{jobname}/analysis/pca')

    # Display configurations
    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Predictions Path: {predictions_path}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")
    print(f"Analysis Range: {analysis_range_name} = {analysis_range}")
    print(f"Align Range: {align_range}")
    print(f"Number of Clusters: {n_pca_clusters}")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    print("***************************************************************\n")

    # Load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # Run PCA analysis
    pca_from_ensemble(jobname, pre_analysis_dict, output_path, align_range, analysis_range, int(n_pca_clusters), widget)


def main():
    """
    Main function to parse arguments and run PCA analysis.
    """

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Runs PCA analysis on an ensemble of structures.")

    parser.add_argument('--config_file', type=str, default='config.json',
                        help="OPTIONAL: Path to load configuration from file (default: config.json)")
    parser.add_argument('--output_path', type=str,
                        help="Path to read predictions and save results to (default: current directory)")
    parser.add_argument('--predictions_path', type=str,
                        help="Path to read predictions from (default: determined by output_path and jobname)")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--starting_residue', type=int,
                        help="Sets the starting residue for reindexing (predictions are usually 1-indexed)")
    parser.add_argument('--align_range', type=str, help="The atom alignment range for PCA calculations (MDAnalysis Syntax)")
    parser.add_argument('--seq_pairs', type=str, help="The sequence pairs used to generate the predictions")
    parser.add_argument('--analysis_range', type=str, help="The atom range for PCA calculations after alignment to align_range")
    parser.add_argument('--analysis_range_name', type=str, help="The name of the atom range (e.g. kinase core, helix 1, etc.)")
    parser.add_argument('--engine', type=str, help="The engine used to generate predictions (e.g. AlphaFold2, OpenFold)")
    parser.add_argument('--n_pca_clusters', type=int, help="Number of PCA clusters to generate")

    # Parse arguments
    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    config.update({k: v for k, v in vars(args).items() if v is not None})

    # Run PCA analysis with the provided configuration
    run_pca_analysis(config)


if __name__ == "__main__":
    main()
