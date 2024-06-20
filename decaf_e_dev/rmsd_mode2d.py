import os
import argparse
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

from decaf_e_dev.ensemble_analysis.analysis_utils import create_directory
from decaf_e_dev.ensemble_analysis.analysis_utils import load_predictions
from decaf_e_dev.ensemble_analysis.analysis_utils import load_config
from decaf_e_dev.ensemble_analysis.analysis_utils import auto_select_2d_references
from decaf_e_dev.ensemble_analysis.twodrmsd import TwodRMSD

## TODO: add sanity checks, rewrite help snipets

def main():
    parser = argparse.ArgumentParser(description="Runs state detection analysis using RMSD vs. two references (2D) "
                                                 "previously identified with rmsd_mode1d.py")

    parser.add_argument('--config_file', type=str, help="Path to the configuration file")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--output_path', type=str, help="Path to save results to")
    parser.add_argument('--mode_results', type=str,
                        help="Path to output .csv file from 1D mode/state detection tool "
                             "(if not provided, will search automatically based on other parameters)")
    parser.add_argument('--seq_pairs', type=str, help="A list of [max_seq, extra_seq] pairs "
                                                      "previously used to generate the predictions")
    parser.add_argument('--predictions_path', type=str,
                        help="Path to read PDB files of predictions, expects format: /jobname_maxseq_extraseq/ "
                             "(if not provided, will search automatically based on other parameters)")
    parser.add_argument('--engine', type=str, choices=['alphafold2', 'openfold'],
                        help="The engine previously used to generate predictions (AlphaFold2 or OpenFold), "
                             "used to find predictions if predictions_path is not supplied")
    parser.add_argument('--starting_residue', type=int,
                        help="Sets the starting residue for reindexing (predictions are usually 1-indexed)")
    parser.add_argument('--align_range', type=str, help="The atom alignment range for RMSF calculations "
                                                        "(MDAnalysis Syntax)")
    parser.add_argument('--analysis_range', type=str, help="The atom range for RMSD calculations "
                                                           "after alignment to --align_range")
    parser.add_argument('--analysis_range_name', type=str, help="The name of the atom range "
                                                                "(e.g. kinase core, helix 1, etc.)")
    parser.add_argument('--ref2d1', type=str, help="Path to the .pdb file defining "
                                                   "the first reference structure "
                                                  "for RMSD calculations "
                                                  "(if not provided, chooses from mode dataset)")
    parser.add_argument('--ref2d2', type=str, help="Path to the .pdb file defining"
                                                   "the second reference structure "
                                                  "for RMSD calculations "
                                                  "(if not provided, chooses from mode datasetn)")
    parser.add_argument('--n_stdevs', type=str, help="Number of standard deviations "
                                                     "to consider when calculating close points to fit curve")
    parser.add_argument('--n_clusters', type=str, help="Number of standard deviations "
                                                     "to consider when calculating close points to fit curve")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    jobname = args.jobname if args.jobname else config.get('jobname')
    output_path = args.output_path if args.output_path else config.get('output_path')
    mode_results = args.mode_results if args.mode_results else config.get('mode_results')
    seq_pairs = args.seq_pairs if args.seq_pairs else config.get('seq_pairs')
    predictions_path = args.predictions_path if args.predictions_path else config.get('predictions_path')
    engine = args.engines if args.engine else config.get('engine')
    align_range = args.align_range if args.align_range else config.get('align_range')
    analysis_range = args.analysis_range if args.analysis_range else config.get('analysis_range')
    analysis_range_name = args.analysis_range_name if args.analysis_range_name else config.get('analysis_range_name')
    ref2d1 = args.ref1 if args.ref2d1 else config.get('ref2d1')
    ref2d2 = args.ref2 if args.ref2d2 else config.get('ref2d2')
    n_stdevs = args.n_stdevs if args.n_stdevs else config.get('n_stdevs')
    n_clusters = args.n_clusters if args.n_clusters else config.get('n_clusters')
    starting_residue = args.starting_residue if args.starting_residue else config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    if not mode_results:
        mode_results = f'{output_path}/{jobname}/analysis/rmsd_1d/{jobname}_rmsd_1d_analysis_results.csv'

    if not ref2d1 and not ref2d2:
        ref2d1, ref2d2 = auto_select_2d_references(mode_results, 'RMSD')

    create_directory(f'{output_path}/{jobname}/analysis/rmsd_2d')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Used Config File? {args.config_file is not None}")
    print(f"Job Name: {jobname}")
    print(f"Output Path: {output_path}")
    print(f"max_seq:extra_seq Pairs: {seq_pairs}")
    print(f"Predictions Path: {predictions_path}")
    print(f"Engine: {engine}")
    print(f"Align Range: {align_range}")
    print(f"Analysis Range: {analysis_range_name} = {analysis_range}")
    print(f"Reference 1: {ref2d1}")
    print(f"Reference 2: {ref2d2}")
    print(f"Number of Standard Devs. to Consider Points Close to fit Line: {n_stdevs}")
    if n_clusters:
        print(f"Number of Clusters: {n_clusters}")
    else:
        print(f"Number of Clusters: Number of Detected 1D RMSD Modes + 1")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    print("***************************************************************\n")

    input_dict = {'jobname': jobname,
                  'output_path': output_path,
                  'seq_pairs': seq_pairs,
                  'analysis_range': analysis_range,
                  'analysis_range_name': analysis_range_name,
                  'align_range': align_range,}

    # load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # rum 2D RMSD analysis
    twod = TwodRMSD(pre_analysis_dict, input_dict, ref2d1, ref2d2)

    # build and save results dataset
    twod.get_2d_rmsd(mode_results, n_stdevs, n_clusters)


if __name__ == "__main__":
    main()
