import os
import argparse

from ensemble_analysis.analysis_utils import create_directory
from ensemble_analysis.analysis_utils import load_predictions
from ensemble_analysis.analysis_utils import load_config
from ensemble_analysis.traj import save_trajs
def main():

    parser = argparse.ArgumentParser(description="Builds trajectory from PDB predictions, "
                                                 "optionally reordering based on a dataset to be sorted")

    parser.add_argument('--config_file', type=str, default='config.json',
                        help="OPTIONAL: Path to load configuration from file (default: config.json)")
    parser.add_argument('--output_path', type=str,
                        help="Path to read AF2 predictions and save results to (default: current directory)")
    parser.add_argument('--predictions_path', type=str,
                        help="Path to read AF2 predictions and save results to (default: current directory)")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--seq_pairs', type=str, help="The job name")
    parser.add_argument('--analysis_range', type=str, help="The job name")
    parser.add_argument('--analysis_range_name', type=str, help="The job name")
    parser.add_argument('--reorder', type=str, help="The job name")
    parser.add_argument('--traj_format', type=str, help="The job name")
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
    analysis_range = args.analysis_range if args.analysis_range else config.get('analysis_range')
    analysis_range_name = args.analysis_range_name if args.analysis_range_name else config.get('analysis_range_name')
    reorder = args.reorder if args.reorder else config.get('reorder')
    traj_format = args.traj_format if args.traj_format else config.get('traj_format')
    engine = args.engines if args.engine else config.get('engine')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    create_directory(f'{output_path}/{jobname}/trajs')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Used Config File? {args.config_file is not None}")
    print(f"Predictions Path: {predictions_path}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")
    print(f"Analysis Range: {analysis_range_name} = {analysis_range}")
    print(f"Reorder by: {reorder}")
    print(f"Traj format: {traj_format}")
    print(f"Engine: {engine}")
    print("***************************************************************\n")

    input_dict = {'jobname': jobname,
                  'output_path': output_path,
                  'seq_pairs': seq_pairs,
                  'analysis_range': analysis_range,
                  'analysis_range_name': analysis_range_name}

    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname)

    save_trajs(pre_analysis_dict, input_dict, reorder, traj_format)


if __name__ == "__main__":
    main()
