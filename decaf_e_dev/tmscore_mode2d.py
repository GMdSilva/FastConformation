import os
import argparse

import warnings
warnings.filterwarnings("ignore")

from ensemble_analysis.analysis_utils import create_directory
from ensemble_analysis.analysis_utils import load_predictions
from ensemble_analysis.analysis_utils import load_config
from ensemble_analysis.analysis_utils import auto_select_2d_references
from ensemble_analysis.twotmscore import TwoTMScore


def main():
    parser = argparse.ArgumentParser(description="Build jackhmmer MSA for a list of sequences.")

    parser.add_argument('--config_file', type=str, default='config.json',
                        help="OPTIONAL: Path to load configuration from file (default: config.json)")
    parser.add_argument('--output_path', type=str,
                        help="Path to read AF2 predictions and save results to (default: current directory)")
    parser.add_argument('--predictions_path', type=str,
                        help="Path to read AF2 predictions and save results to (default: current directory)")
    parser.add_argument('--mode_results', type=str,
                        help="Path to read AF2 predictions and save results to (default: current directory)")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--seq_pairs', type=str, help="The job name")
    parser.add_argument('--starting_residue', type=int,
                        help="Sets the starting residue for reindexing (predictions are usually 1-indexed)")
    parser.add_argument('--slice_predictions', type=str, help="The job name")
    parser.add_argument('--engine', type=str, help="The job name")
    parser.add_argument('--ref2d1', type=str, help="The job name")
    parser.add_argument('--ref2d2', type=str, help="The job name")
    parser.add_argument('--n_stdevs', type=str, help="The job name")
    parser.add_argument('--n_clusters', type=str, help="Number of standard deviations "
                                                     "to consider when calculating close points to fit curve")
    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    output_path = args.output_path if args.output_path else config.get('output_path')
    predictions_path = args.predictions_path if args.predictions_path else config.get('predictions_path')
    mode_results = args.mode_results if args.mode_results else config.get('mode_results')
    seq_pairs = args.seq_pairs if args.seq_pairs else config.get('seq_pairs')
    jobname = args.jobname if args.jobname else config.get('jobname')
    ref2d1 = args.ref1 if args.ref2d1 else config.get('ref2d1')
    ref2d2 = args.ref2 if args.ref2d2 else config.get('ref2d2')
    starting_residue = args.starting_residue if args.starting_residue else config.get('starting_residue')
    slice_predictions = args.slice_predictions if args.slice_predictions else config.get('slice_predictions')
    engine = args.engines if args.engine else config.get('engine')
    n_stdevs = args.n_stdevs if args.n_stdevs else config.get('n_stdevs')
    n_clusters = args.n_clusters if args.n_clusters else config.get('n_clusters')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    if not mode_results:
        mode_results = f'{output_path}/{jobname}/analysis/tmscore_1d/{jobname}_tmscore_1d_analysis_results.csv'

    if not ref2d1 and not ref2d2:
        ref2d1, ref2d2 = auto_select_2d_references(mode_results, 'tmscore')

    create_directory(f'{output_path}/{jobname}/analysis/tmscore_2d')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Used Config File? {args.config_file is not None}")
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

    input_dict = {'jobname': jobname,
                  'output_path': output_path,
                  'seq_pairs': seq_pairs,
                  'predictions_path': predictions_path}

    # loads predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # runs 2D TM-Score analysis
    twod = TwoTMScore(pre_analysis_dict, input_dict, ref2d1, ref2d2, slice_predictions)

    # builds results dataset and saves to disk
    twod.get_2d_tmscore(mode_results, n_stdevs, n_clusters)


if __name__ == "__main__":
    main()
