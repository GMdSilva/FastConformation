import os
import argparse

import warnings
warnings.filterwarnings("ignore")

from ensemble_analysis.analysis_utils import create_directory
from ensemble_analysis.analysis_utils import load_predictions
from ensemble_analysis.analysis_utils import load_config
from ensemble_analysis.rmsd import rmsd_mode_analysis, build_dataset_rmsd_modes

## TODO: add sanity checks, rewrite help snipets

def main():
    parser = argparse.ArgumentParser(description="Runs state detection analysis using RMSD vs. a single reference (1D)")

    parser.add_argument('--config_file', type=str, help="Path to the configuration file")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--output_path', type=str, help="Path to save results to")
    parser.add_argument('--seq_pairs', type=str, help="A list of [max_seq, extra_seq] pairs "
                                                      "previously used to generate the predictions")
    parser.add_argument('--predictions_path', type=str,
                        help="Path to read PDB files of predictions, expects format: /jobname_maxseq_extraseq/ "
                             "(if not provided, will search automatically based on other parameters)")
    parser.add_argument('--engine', type=str, choices=['alphafold2', 'openfold'],
                        help="The engine previously used to generate predictions (AlphaFold2 or OpenFold), "
                             "used to find predictions if predictions_path is not supplied")
    parser.add_argument('--align_range', type=str, help="The atom alignment range for RMSF calculations "
                                                        "(MDAnalysis Syntax)")
    parser.add_argument('--analysis_range', type=str, help="The atom range for RMSD calculations "
                                                           "after alignment to --align_range")
    parser.add_argument('--analysis_range_name', type=str, help="The name of the atom range "
                                                                "(e.g. kinase core, helix 1, etc.)")
    parser.add_argument('--ref1d', type=str, help="Path to the .pdb file defining the reference structure "
                                                  "for RMSD calculations "
                                                  "(if not provided, picks top ranked prediction)")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    jobname = args.jobname if args.jobname else config.get('jobname')
    output_path = args.output_path if args.output_path else config.get('output_path')
    seq_pairs = args.seq_pairs if args.seq_pairs else config.get('seq_pairs')
    predictions_path = args.predictions_path if args.predictions_path else config.get('predictions_path')
    engine = args.engines if args.engine else config.get('engine')
    align_range = args.align_range if args.align_range else config.get('align_range')
    analysis_range = args.analysis_range if args.analysis_range else config.get('analysis_range')
    analysis_range_name = args.analysis_range_name if args.analysis_range_name else config.get('analysis_range_name')
    ref1d = args.ref1d if args.ref1d else config.get('ref1d')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    create_directory(f'{output_path}/{jobname}/analysis/rmsd_1d')
    create_directory(f'{output_path}/{jobname}/analysis/representative_structures/rmsd_1d')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Used Config File? {args.config_file is not None}")
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
    print("***************************************************************\n")

    input_dict = {'jobname': jobname,
                  'output_path': output_path,
                  'seq_pairs': seq_pairs,
                  'analysis_range': analysis_range,
                  'analysis_range_name': analysis_range_name,
                  'align_range': align_range}

    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname)
    rmsd_mode_analysis_dict = rmsd_mode_analysis(pre_analysis_dict, input_dict, ref1d)
    build_dataset_rmsd_modes(rmsd_mode_analysis_dict, input_dict)


if __name__ == "__main__":
    main()
