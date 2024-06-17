import subprocess
import os
import argparse
import json
import shutil

from msa_generation.msa_utils import read_fasta, save_dict_to_fasta

def validate_inputs(fasta_path, output_path, jobname):
    # Check if fasta_path is a valid file
    if not os.path.isfile(fasta_path):
        raise ValueError(f"MSA path '{fasta_path}' is not a valid file.")

    if not os.path.isdir(output_path):
        raise ValueError(f"Output path '{output_path}' is not a valid directory.")

    # Check if jobname is a string
    if not isinstance(jobname, str):
        raise ValueError(f"Jobname '{jobname}' is not a valid string.")


def get_mmseqs_msa(sequence_file_path, sequence_name, output_path, jobname, env):
    complete_output_path = f'{output_path}/{jobname}'

    command = (f"colabfold_batch {sequence_file_path} {complete_output_path} --msa-only")

    print(f'Getting mmseqs2 MSA for {jobname}')

    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)

    # Read line by line as they are output
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()

    copy_msa_and_clean(f"{complete_output_path}/{sequence_name}.a3m", output_path)


def load_config(config_file):
    # Default configuration values
    default_config = {
        'fasta_path': None,
        'output_path': './',
        'jobname': 'prediction_run',
        'seq_pairs': [[8, 16], [64, 128], [512, 1024]],
        'seeds': 10,
        'save_all': False,
        'platform': 'cpu'
    }

    try:
        with open(config_file, 'r') as file:
            config = json.load(file)
        return config
    except FileNotFoundError:
        print(f"Configuration file {config_file} not found. Using default values.")
        return default_config
    except json.JSONDecodeError:
        print(f"Error decoding JSON from the configuration file {config_file}. Using default values.")
        return default_config


def save_config(config, file_path):
    with open(file_path, 'w') as file:
        json.dump(config, file, indent=4)


def remove_first_line(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        file.writelines(lines[1:])


def remove_duplicate_sequences(input_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    seq_dict = {}
    current_name = None
    counter = 0
    for line in lines:
        if line.startswith('>'):
            current_name = line.strip()
        else:
            seq = line.strip()
            if seq not in seq_dict:
                seq_dict[seq] = current_name
            else:
                counter += 1

    with open(input_file, 'w') as f:
        for seq, name in seq_dict.items():
            f.write(f"{name}\n{seq}\n")


def copy_msa_and_clean(src_file_path, dest_dir_path):
    remove_first_line(src_file_path)
    remove_duplicate_sequences(src_file_path)
    # Ensure the destination directory exists
    os.makedirs(dest_dir_path, exist_ok=True)

    # Copy the file to the destination directory
    shutil.copy(src_file_path, dest_dir_path)

    # Get the directory of the source file
    src_dir_path = os.path.dirname(src_file_path)

    # Delete the source directory and all its contents
    shutil.rmtree(src_dir_path)


def main():
    parser = argparse.ArgumentParser(description="Run multiple ensemble predictions.")

    parser.add_argument('--config_file', type=str, default='config.json',
                        help="Path to the configuration file (default: config.json)")
    parser.add_argument('--fasta_path', type=str,
                        help="Path to the .fasta file containing target sequences to build the MSA for")
    parser.add_argument('--output_path', type=str, help="Path to save results to"
                                                        "(default: current directory)")
    parser.add_argument('--jobname', type=str, help="The job name")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    fasta_path = args.fasta_path if args.fasta_path else config.get('fasta_path')
    output_path = args.output_path if args.output_path else config.get('output_path')
    jobname = args.jobname if args.jobname else config.get('jobname')
    platform = 'cpu'

    os.environ['JAX_PLATFORMS'] = platform

    env = os.environ.copy()
    env["PATH"] += os.pathsep + "localcolabfold/colabfold-conda/bin"
    env["PATH"] += os.pathsep + "/home/gabriel/localcolabfold/colabfold-conda/bin"  ## TODO: remove path

    # Print out the configurations for debugging
    print("Configurations:")
    print(f"Used Config File? {args.config_file is not None}")
    print(f"Sequence file path: {fasta_path}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")

    sequences = read_fasta(fasta_path)
    save_dict_to_fasta(sequences, output_path)

    for sequence in sequences:
        single_fasta_path = f'{output_path}/{sequence}.fasta'
        get_mmseqs_msa(single_fasta_path, sequence, output_path, jobname, env)

    save_config(config, config_file)


if __name__ == "__main__":
    main()
