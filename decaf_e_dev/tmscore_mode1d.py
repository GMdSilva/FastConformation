import os
import argparse

import warnings
warnings.filterwarnings("ignore")

from decaf_e_dev.ensemble_analysis.analysis_utils import create_directory, load_predictions
from decaf_e_dev.ensemble_analysis.analysis_utils import load_config
from decaf_e_dev.ensemble_analysis.tmscore import build_dataset_tmscore_modes, tmscore_mode_analysis


def main():
    parser = argparse.ArgumentParser(description="Build jackhmmer MSA for a list of sequences.")

    parser.add_argument('--config_file', type=str, default='config.json',
                        help="OPTIONAL: Path to load configuration from file (default: config.json)")
    parser.add_argument('--output_path', type=str,
                        help="Path to read AF2 predictions and save results to (default: current directory)")
    parser.add_argument('--predictions_path', type=str,
                        help="Path to read AF2 predictions and save results to (default: current directory)")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--seq_pairs', type=str, help="The job name")
    parser.add_argument('--starting_residue', type=int,
                        help="Sets the starting residue for reindexing (predictions are usually 1-indexed)")
    parser.add_argument('--slice_predictions', type=str, help="The job name")
    parser.add_argument('--ref1', type=str, help="The job name")
    parser.add_argument('--engine', type=str, help="The job name")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    output_path = args.output_path if args.output_path else config.get('output_path')
    predictions_path = args.predictions_path if args.predictions_path else config.get('predictions_path')
    seq_pairs = args.seq_pairs if args.seq_pairs else config.get('seq_pairs')
    jobname = args.jobname if args.jobname else config.get('jobname')
    slice_predictions = args.slice_predictions if args.slice_predictions else config.get('slice_predictions')
    ref1 = args.ref1 if args.ref1 else config.get('ref1')
    engine = args.engines if args.engine else config.get('engine')
    starting_residue = args.starting_residue if args.starting_residue else config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    create_directory(f'{output_path}/{jobname}/analysis/tmscore_1d')
    create_directory(f'{output_path}/{jobname}/analysis/representative_structures/tmscore_1d')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Used Config File? {args.config_file is not None}")
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

    input_dict = {'jobname': jobname,
                  'predictions_path': predictions_path,
                  'output_path': output_path,
                  'seq_pairs': seq_pairs}

    # load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # runs 1D TM-Score analysis
    tmscore_mode_analysis_dict = tmscore_mode_analysis(pre_analysis_dict, input_dict, ref1, slice_predictions)

    # builds results dataset and saves to disk
    build_dataset_tmscore_modes(tmscore_mode_analysis_dict, input_dict)


if __name__ == "__main__":
    main()
