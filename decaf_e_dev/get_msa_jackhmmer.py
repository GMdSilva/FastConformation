import os
import argparse
import json
import pickle

from msa_generation import get_msa_jackhmmer
from msa_generation.msa_utils import create_ram_disk, read_fasta


def load_config(file_path):
    # Default configuration values
    default_config = {
        'sequences_file': None,
        'output_path': './',
        'jobname': "prediction_run",
        'homooligomers': 1,
        'tmp_dir': "./tmp",
        'use_ramdisk': False
    }

    try:
        with open(file_path, 'r') as file:
            config = json.load(file)
        return config
    except FileNotFoundError:
        print(f"Configuration file {file_path} not found. Using default values.")
        return default_config
    except json.JSONDecodeError:
        print(f"Error decoding JSON from the configuration file {file_path}. Using default values.")
        return default_config


def save_config(config, file_path):
    with open(file_path, 'w') as file:
        json.dump(config, file, indent=4)


def reformat_sequences(input_msa):
    formatted_sequences = []
    for idx, sequence in enumerate(input_msa[0]):
        formatted_sequence = f">sequence_{idx + 1}\n{sequence}\n"
        formatted_sequences.append(formatted_sequence)
    return formatted_sequences


def convert_msa(filename):
    with open(filename, 'rb') as f:
        msa_data = pickle.load(f)

    msa_list = msa_data['msas']
    converted_msa = reformat_sequences(msa_list)
    return converted_msa


def save_formatted_sequences_to_file(formatted_sequences, output_file):
    with open(f'{output_file}', 'w') as f:
        for sequence in formatted_sequences:
            f.write(sequence)


def prepare_os(tmp_dir, use_ramdisk):
    os.makedirs(tmp_dir, exist_ok=True)
    if use_ramdisk:
        create_ram_disk()


def build_msa(sequence, jobname, complete_output_dir, homooligomer, tmp_dir, use_ramdisk):
    prepped_msa = get_msa_jackhmmer.prep_inputs(sequence, jobname, homooligomer,  output_dir=complete_output_dir)

    get_msa_jackhmmer.prep_msa(prepped_msa, msa_method='jackhmmer', add_custom_msa=False, msa_format="fas",
                     pair_mode="unpaired", pair_cov=50, pair_qid=20, TMP_DIR=tmp_dir, use_ramdisk=use_ramdisk)


def main():
    parser = argparse.ArgumentParser(description="Build jackhmmer MSA for a list of sequences.")

    parser.add_argument('--config_file', type=str, default='config.json',
                        help="OPTIONAL: Path to load configuration from file (default: config.json)")
    parser.add_argument('--sequences_file', type=str,
                        help="Path to a .fasta file containing target sequences for MSA building.")
    parser.add_argument('--output_path', type=str,
                        help="Path to save results to (default: current directory)")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--homooligomers', type=int,
                        help="Number of copies of the protein, defaults to 1 (monomer)")
    parser.add_argument('--use_ramdisk', type=bool, help="If you have root access,"
                                                         "mounts a ramdisk for dramatically speeding up jackhmmer.")
    parser.add_argument('--tmp_dir', type=str, help="Temporary directory to handle IO (default: ./tmp)")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config_msa.json'
    config = load_config(config_file)

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
