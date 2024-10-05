import os
import warnings
import argparse
from fast_conformation.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config
from fast_conformation.ensemble_analysis.rmsd import rmsd_mode_analysis, build_dataset_rmsd_modes

warnings.filterwarnings("ignore")


def run_rmsd_analysis(config, plot_widget=None):
    """
    Run RMSD analysis based on the provided configuration.

    Parameters:
    config (dict): Configuration dictionary containing parameters for the analysis.
    plot_widget (object, optional): Plot widget for displaying results (default is None).

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
    analysis_range = config.get('analysis_range')
    analysis_range_name = config.get('analysis_range_name')
    ref1d = config.get('ref1d') if config.get('ref1d') else None
    starting_residue = config.get('starting_residue')

    # Check if the output path is a valid directory
    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    # Set default predictions path if not provided
    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    # Create necessary directories
    create_directory(f'{output_path}/{jobname}/analysis/rmsd_1d')
    create_directory(f'{output_path}/{jobname}/analysis/representative_structures/rmsd_1d')

    # Display configurations
    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Output Path: {output_path}")
    print(f"max_seq:extra_seq Pairs: {seq_pairs}")
    print(f"Predictions Path: {predictions_path}")
    print(f"Engine: {engine}")
    print(f"Job Name: {jobname}")
    print(f"Align Range: {align_range}")
    print(f"Analysis Range: {analysis_range_name} = {analysis_range}")
    if ref1d:
        print(f"Reference Structure: {ref1d}")
    else:
        print(f"Reference Structure: Top Ranked Prediction")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    print("***************************************************************\n")

    # Prepare input dictionary
    input_dict = {
        'jobname': jobname,
        'output_path': output_path,
        'seq_pairs': seq_pairs,
        'analysis_range': analysis_range,
        'analysis_range_name': analysis_range_name,
        'align_range': align_range
    }

    # Load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # Run 1D RMSD analysis
    rmsd_mode_analysis_dict = rmsd_mode_analysis(pre_analysis_dict, input_dict, ref1d, plot_widget)

    # Build and save results dataset
    build_dataset_rmsd_modes(rmsd_mode_analysis_dict, input_dict)


def main():
    """
    Main function to parse arguments and run RMSD analysis.
    """

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Runs state detection analysis using RMSD vs. a single reference (1D)")

    parser.add_argument('--config_file', type=str, help="Path to the configuration file")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--output_path', type=str, help="Path to save results to")
    parser.add_argument('--seq_pairs', type=str, help="A list of [max_seq, extra_seq] pairs previously used to generate the predictions")
    parser.add_argument('--predictions_path', type=str, help="Path to read PDB files of predictions, expects format: /jobname_maxseq_extraseq/ (if not provided, will search automatically based on other parameters)")
    parser.add_argument('--engine', type=str, choices=['alphafold2', 'openfold'], help="The engine previously used to generate predictions (AlphaFold2 or OpenFold), used to find predictions if predictions_path is not supplied")
    parser.add_argument('--starting_residue', type=int, help="Sets the starting residue for reindexing (predictions are usually 1-indexed)")
    parser.add_argument('--align_range', type=str, help="The atom alignment range for RMSF calculations (MDAnalysis Syntax)")
    parser.add_argument('--analysis_range', type=str, help="The atom range for RMSD calculations after alignment to --align_range")
    parser.add_argument('--analysis_range_name', type=str, help="The name of the atom range (e.g. kinase core, helix 1, etc.)")
    parser.add_argument('--ref1d', type=str, help="Path to the .pdb file defining the reference structure for RMSD calculations (if not provided, picks top ranked prediction)")

    # Parse arguments
    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    config.update({k: v for k, v in vars(args).items() if v is not None})

    # Run RMSD analysis with the provided configuration
    run_rmsd_analysis(config)


if __name__ == "__main__":
    main()
