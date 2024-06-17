from ensemble_analysis.analysis_utils import load_config
from ensemble_analysis.analysis_utils import create_directory
from ensemble_analysis.analysis_utils import load_predictions, load_predictions_json
from ensemble_analysis.rmsf import calculate_rmsf_multiple, calculate_rmsf_and_call_peaks, build_dataset_rmsf_peaks
#from scripts.pca_analysis_script import pca_from_ensemble_movie

import os
import argparse
import json
import pickle

def main():
    parser = argparse.ArgumentParser(description="Build jackhmmer MSA for a list of sequences.")

    parser.add_argument('--config_file', type=str, default='config.json',
                        help="OPTIONAL: Path to load configuration from file (default: config.json)")
    parser.add_argument('--output_path', type=str,
                        help="Path to save results to (default: current directory)")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--align_range', type=str, help="The job name")
    parser.add_argument('--seq_pairs', type=str, help="The job name")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config_msa.json'
    config = load_config(config_file)

    config = load_config('example_configs/config.json')

    # Override config with command line arguments if provided
    output_path = config.get('output_path')
    seq_pairs = config.get('seq_pairs')
    jobname = config.get('jobname')
    align_range = config.get('align_range')
    #
    create_directory(f'{output_path}/representative_structures')
    create_directory(f'{output_path}/plots')
    create_directory(f'{output_path}/datasets')

    input_dict = {'jobname': jobname,
                  'output_path': output_path,
                  'seq_pairs': seq_pairs,
                  'rmsd_range': analysis_range,
                  'rmsd_range_name': analysis_range_name,
                  'align_range': align_range, }
    pre_analysis_dict = load_predictions(output_path, seq_pairs, jobname)

    calculate_rmsf_multiple(jobname, pre_analysis_dict, align_range, output_path)
    rmsf_peak_calling_dict = calculate_rmsf_and_call_peaks(jobname, pre_analysis_dict, align_range, output_path)
    build_dataset_rmsf_peaks(input_dict['jobname'], rmsf_peak_calling_dict, output_path)


    # Override config with command line arguments if provided
    sequences_file = args.sequences_file if args.sequences_file else config.get('sequences_file')
    output_path = args.output_path if args.output_path else config.get('output_path')
    jobname = args.jobname if args.jobname else config.get('jobname')
    homooligomers = args.homooligomers if args.homooligomers else config.get('homooligomers')
    tmp_dir = args.tmp_dir if args.tmp_dir else config.get('tmp_dir')
    use_ramdisk = args.use_ramdisk if args.use_ramdisk else config.get('use_ramdisk')

    if sequences_file is None:
        raise ValueError("Sequences file must be provided via --sequences_file or in the config file.")
    if not os.path.exists(sequences_file):
        raise FileNotFoundError(f"Sequences file {sequences_file} not found.")
    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory.")
    if not isinstance(homooligomers, int) or homooligomers <= 0:
        raise ValueError("Homooligomers must be a positive integer.")

    print("Configurations:")
    print(f"Used Config File? {args.config_file is not None}")
    print(f"Sequences Path: {sequences_file}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")
    print(f"Homooligomer: {homooligomers}")
    print(f"Temporary Directory: {tmp_dir}")
    print(f"Use ramdisk? {use_ramdisk}")

    prepare_os(tmp_dir, use_ramdisk)

    sequences = read_fasta(sequences_file)

    for sequence in sequences:
        sequence_name = sequence
        sequence_string = sequences[sequence]

        print(f"Building jackhmmer MSA for {sequence_name}, with sequence\n{sequence_string}")

        complete_output_dir = f"{output_path}/{sequence_name}"

        build_msa(sequence_string, jobname, complete_output_dir, homooligomers, tmp_dir, use_ramdisk)
        converted_msa = convert_msa(f"{complete_output_dir}/msa.pickle")
        save_formatted_sequences_to_file(converted_msa, f"{complete_output_dir}/{jobname}_{sequence_name}.a3m")

    save_config(config, config_file)


if __name__ == "__main__":
    main()

