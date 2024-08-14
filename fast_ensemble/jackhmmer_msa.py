import os
import argparse
import json
import pickle

from fast_ensemble.msa_generation import get_msa_jackhmmer
from fast_ensemble.msa_generation.msa_utils import create_ram_disk, read_fasta, create_directory, save_dict_to_fasta


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


def prepare_os():
    create_ram_disk()


def build_msa(sequence, jobname, complete_output_dir, homooligomer, tmp_dir, use_ramdisk):
    prepped_msa = get_msa_jackhmmer.prep_inputs(sequence, jobname, homooligomer,  output_dir=complete_output_dir)

    get_msa_jackhmmer.prep_msa(prepped_msa, msa_method='jackhmmer', add_custom_msa=False, msa_format="fas",
                     pair_mode="unpaired", pair_cov=50, pair_qid=20, TMP_DIR=tmp_dir, use_ramdisk=use_ramdisk)


def main():
    parser = argparse.ArgumentParser(description="Assemble MSA for target sequence with jackhmmer")
    parser.add_argument('--config_file', type=str, help="Path to the configuration file")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--sequence_path', type=str,
                        help="Path to a .fasta file containing the target sequence for MSA building")
    parser.add_argument('--output_path', type=str, help="Path to save results to")
    parser.add_argument('--homooligomers', type=int, help="Number of copies of the protein")
    parser.add_argument('--use_ramdisk', type=bool, help="Mounts a ramdisk for dramatically "
                                                         "speeding up jackhmmer search (requires root access)")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    jobname = args.jobname if args.jobname else config.get('jobname')
    sequence_path = args.sequence_path if args.sequence_path else config.get('sequence_path')
    output_path = args.output_path if args.output_path else config.get('output_path')
    homooligomers = args.homooligomers if args.homooligomers else config.get('homooligomers')
    use_ramdisk = args.use_ramdisk if args.use_ramdisk else config.get('use_ramdisk')

    if sequence_path is None:
        raise ValueError("Sequences file must be provided via --sequences_file or in the config file.")
    if not os.path.exists(sequence_path):
        raise FileNotFoundError(f"Sequences file {sequence_path} not found.")
    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory.")
    if not isinstance(homooligomers, int) or homooligomers <= 0:
        raise ValueError("Homooligomers must be a positive integer.")

    create_directory(f'{output_path}/{jobname}/msas/jackhmmer')
    create_directory(f'{output_path}/{jobname}/target_seq/')

    sequence_dict = read_fasta(sequence_path)
    save_dict_to_fasta(sequence_dict, output_path, jobname)

    sequence_name = list(sequence_dict.keys())[0]
    sequence_string = list(sequence_dict.values())[0]

    print(f"\nRunning jackhmmer MSA building for target sequence at {sequence_path}\n")

    print(f"{sequence_name}: \n{sequence_string}\n")

    print("Configurations:")
    print("***************************************************************")
    print(f"Used Config File? {args.config_file is not None}")
    print(f"Job Name: {jobname}")
    print(f"Sequence File Path: {sequence_path}")
    print(f"Output Path: {output_path}")
    print(f"Homooligomers: {homooligomers}")
    print(f"Use ramdisk? {use_ramdisk}")
    print("***************************************************************\n")

    # optionally create ramdisk
    if use_ramdisk:
        prepare_os()

    # sets output directory for MSA
    complete_output_dir = f"{output_path}/{jobname}/msas/jackhmmer/"

    # builds jackhmmer MSA from target sequence
    build_msa(sequence_string, jobname, complete_output_dir, homooligomers, tmp_dir='tmp', use_ramdisk=use_ramdisk)

    # reformats MSA to something colabfold_batch can use
    converted_msa = convert_msa(f"{complete_output_dir}/msa.pickle")
    save_formatted_sequences_to_file(converted_msa, f"{complete_output_dir}/{jobname}.a3m")

    print(f'\nSaved {jobname} jackhmmer MSA to {complete_output_dir}/{jobname}.a3m\n')


if __name__ == "__main__":
    main()
