import subprocess
import os
import argparse
import json
import shutil

from fast_ensemble.msa_generation.msa_utils import read_fasta, save_dict_to_fasta, create_directory

def validate_inputs(fasta_path, output_path, jobname):
    # Check if fasta_path is a valid file
    if not os.path.isfile(fasta_path):
        raise ValueError(f"MSA path '{fasta_path}' is not a valid file.")

    if not os.path.isdir(output_path):
        raise ValueError(f"Output path '{output_path}' is not a valid directory.")

    # Check if jobname is a string
    if not isinstance(jobname, str):
        raise ValueError(f"Jobname '{jobname}' is not a valid string.")


def get_mmseqs_msa(sequence_file_path, output_path, jobname, env):
    complete_output_path = f'{output_path}/{jobname}/msas/mmseqs2'

    command = (f"colabfold_batch {sequence_file_path} {complete_output_path}/temp_msa --msa-only")

    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)

    # Read line by line as they are output
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()

    shutil.copy(f"{complete_output_path}/temp_msa/{jobname}.pickle", complete_output_path)
    shutil.copy(f"{complete_output_path}/temp_msa/{jobname}_coverage.png", complete_output_path)
    copy_msa_and_clean(f"{complete_output_path}/temp_msa/{jobname}.a3m", complete_output_path)

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

    # Copy the file to the destination directory
    shutil.copy(src_file_path, dest_dir_path)

    # Get the directory of the source file
    src_dir_path = os.path.dirname(src_file_path)

    # Delete the source directory and all its contents
    shutil.rmtree(src_dir_path)


def main():
    parser = argparse.ArgumentParser(description="Assemble MSA for target sequence with mmseqs2")
    parser.add_argument('--config_file', type=str, help="Path to the configuration file")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--sequence_path', type=str,
                        help="Path to a .fasta file containing the target sequence for MSA building")
    parser.add_argument('--output_path', type=str,
                        help="Path to save results to")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    jobname = args.jobname if args.jobname else config.get('jobname')
    sequence_path = args.sequence_path if args.sequence_path else config.get('sequence_path')
    output_path = args.output_path if args.output_path else config.get('output_path')
    platform = 'cpu'

    os.environ['JAX_PLATFORMS'] = platform

    env = os.environ.copy()
    env["PATH"] += os.pathsep + "localcolabfold/colabfold-conda/bin"
    env["PATH"] += os.pathsep + "/home/gabriel/localcolabfold/colabfold-conda/bin"  ## TODO: remove path

    create_directory(f'{output_path}/{jobname}/msas/mmseqs2')
    create_directory(f'{output_path}/{jobname}/target_seq/')

    sequence_dict = read_fasta(sequence_path)
    save_dict_to_fasta(sequence_dict, output_path, jobname)

    sequence_name = list(sequence_dict.keys())[0]
    sequence_string = list(sequence_dict.values())[0]

    print(f"\nRunning mmseqs2 MSA building for target sequence at {sequence_path}\n")

    print(f"{sequence_name}: \n{sequence_string}\n")

    print("Configurations:")
    print("***************************************************************")
    print(f"Used Config File? {args.config_file is not None}")
    print(f"Job Name: {jobname}")
    print(f"Sequence File path: {sequence_path}")
    print(f"Output Path: {output_path}")
    print("***************************************************************")

    # loads target sequence
    single_sequence_path = f'{output_path}/{jobname}/target_seq/{jobname}.fasta'

    # gets mmseqs2 MSA
    get_mmseqs_msa(single_sequence_path, output_path, jobname, env)

    print(f'\nSaved {jobname} mmseqs2 MSA to {output_path}/{jobname}/msas/mmseqs2/{jobname}.a3m\n')

if __name__ == "__main__":
    main()
