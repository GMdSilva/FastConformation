import os
import argparse

import warnings
warnings.filterwarnings("ignore")

from fast_ensemble.ensemble_analysis.analysis_utils import (load_config,
                                              create_directory,
                                              load_predictions,
                                              load_predictions_json)
from fast_ensemble.ensemble_analysis.rmsf import (calculate_rmsf_multiple,
                                    calculate_rmsf_and_call_peaks,
                                    build_dataset_rmsf_peaks,
                                    plot_plddt_line,
                                    plot_plddt_rmsf_corr)

def main():
    parser = argparse.ArgumentParser(description="Runs RMSF/pLDDT Analysis for a set of predictions")

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
    parser.add_argument('--starting_residue', type=int,
                        help="Sets the starting residue for reindexing (predictions are usually 1-indexed)")
    parser.add_argument('--align_range', type=str, help="The atom alignment range for RMSF calculations "
                                                        "(MDAnalysis Syntax)")
    parser.add_argument('--detect_mobile', type=bool, help="Pass True to detect mobile residue ranges")
    parser.add_argument('--peak_width', type=int, help="Sets the RMSF peak width threshold "
                                                       "for mobile residue range detection")
    parser.add_argument('--peak_prominence', type=int, help="Sets the RMSF "
                                                                       "peak prominence threshold "
                                                            "for mobile residue range detection")
    parser.add_argument('--peak_height', type=int, help="Sets the RMSF peak height threshold "
                                                            "for mobile residue range detection")


    # TODO ALLOW USER TO OPTIMIZE PEAK_WIDTH, PROMINENCE, AND HEIGHT ON THE FLY

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
    detect_mobile = args.detect_mobile if args.detect_mobile else config.get('detect_mobile')
    peak_width = args.peak_width if args.peak_width else config.get('peak_width')
    peak_prominence = args.peak_prominence if args.peak_prominence else config.get('peak_prominence')
    peak_height = args.peak_height if args.peak_height else config.get('peak_height')
    starting_residue = args.starting_residue if args.starting_residue else config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    create_directory(f'{output_path}/{jobname}/analysis/rmsf_plddt')

    if detect_mobile:
        create_directory(f'{output_path}/{jobname}/analysis/mobile_detection')

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Used Config File? {args.config_file is not None}")
    print(f"Job Name: {jobname}")
    print(f"Output Path: {output_path}")
    print(f"max_seq:extra_seq Pairs: {seq_pairs}")
    print(f"Predictions Path: {predictions_path}")
    print(f"Engine: {engine}")
    print(f"Align Range: {align_range}")
    print(f"Detect Mobile Segments? {detect_mobile}")
    if detect_mobile:
        print(f"Peak Width: {peak_width}")
        print(f"Peak Prominence: {peak_prominence}")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    print("***************************************************************\n")

    # load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # runs RMSF analysis for all predictions
    calculate_rmsf_multiple(jobname, pre_analysis_dict, align_range, output_path)

    # gets pLDDTs for all predictions
    plddt_dict = load_predictions_json(predictions_path, seq_pairs, jobname)
    plot_plddt_line(jobname, plddt_dict, output_path, starting_residue)

    # runs pLDDT/RMSF correlation for all predictions
    plot_plddt_rmsf_corr(jobname, pre_analysis_dict, plddt_dict, output_path)

    # runs mobile residue range detection with RMSF data and save dataset to disk
    if detect_mobile:
        print(f'\nAutomatically detecting mobile segments')
        rmsf_peak_calling_dict = calculate_rmsf_and_call_peaks(jobname,
                                                               pre_analysis_dict,
                                                               align_range,
                                                               output_path,
                                                               peak_width,
                                                               peak_prominence,
                                                               peak_height)

        build_dataset_rmsf_peaks(jobname, rmsf_peak_calling_dict, output_path, engine)

if __name__ == "__main__":
    main()

