import os
import argparse

import warnings
warnings.filterwarnings("ignore")

from ensemble_analysis.analysis_utils import create_directory
from ensemble_analysis.analysis_utils import load_predictions
from ensemble_analysis.analysis_utils import load_config
from ensemble_analysis.pca import pca_from_ensemble


def main():
    parser = argparse.ArgumentParser(description="Build jackhmmer MSA for a list of sequences.")

    parser.add_argument('--config_file', type=str, default='config.json',
                        help="OPTIONAL: Path to load configuration from file (default: config.json)")
    parser.add_argument('--output_path', type=str,
                        help="Path to read AF2 predictions and save results to (default: current directory)")
    parser.add_argument('--predictions_path', type=str,
                        help="Path to read AF2 predictions and save results to (default: current directory)")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--starting_residue', type=int,
                        help="Sets the starting residue for reindexing (predictions are usually 1-indexed)")
    parser.add_argument('--align_range', type=str, help="The job name")
    parser.add_argument('--seq_pairs', type=str, help="The job name")
    parser.add_argument('--analysis_range', type=str, help="The job name")
    parser.add_argument('--analysis_range_name', type=str, help="The job name")
    parser.add_argument('--engine', type=str, help="The job name")
    parser.add_argument('--n_pca_clusters', type=str, help="The job name")


    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    predictions_path = args.predictions_path if args.predictions_path else config.get('predictions_path')
    output_path = args.output_path if args.output_path else config.get('output_path')
    seq_pairs = args.seq_pairs if args.seq_pairs else config.get('seq_pairs')
    jobname = args.jobname if args.jobname else config.get('jobname')
    align_range = args.align_range if args.align_range else config.get('align_range')
    analysis_range = args.analysis_range if args.analysis_range else config.get('analysis_range')
    analysis_range_name = args.analysis_range_name if args.analysis_range_name else config.get('analysis_range_name')
    engine = args.engines if args.engine else config.get('engine')
    n_pca_clusters = args.n_pca_clusters if args.n_pca_clusters else config.get('n_pca_clusters')
    starting_residue = args.starting_residue if args.starting_residue else config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    create_directory(f'{output_path}/{jobname}/analysis/pca')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Used Config File? {args.config_file is not None}")
    print(f"Predictions Path: {predictions_path}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")
    print(f"Analysis Range: {analysis_range_name} = {analysis_range}")
    print(f"Align Range: {align_range}")
    print(f"Number of Clusters: {n_pca_clusters}")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    print("***************************************************************\n")

    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)
    pca_from_ensemble(jobname, pre_analysis_dict, output_path, align_range, analysis_range, int(n_pca_clusters))


if __name__ == "__main__":
    main()
