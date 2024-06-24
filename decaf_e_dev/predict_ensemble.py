import subprocess
import os
import argparse
import json

from decaf_e_dev.msa_generation.msa_utils import create_directory


def validate_inputs(msa_path, output_path, jobname, seq_pairs, seeds):
    # Check if msa_path is a valid file
    if not os.path.isfile(msa_path):
        raise ValueError(f"MSA path '{msa_path}' is not a valid file.")

    if not os.path.isdir(output_path):
        raise ValueError(f"Output path '{output_path}' is not a valid directory.")

    # Check if jobname is a string
    if not isinstance(jobname, str):
        raise ValueError(f"Jobname '{jobname}' is not a valid string.")

    # Check if seq_pairs is a list of pairs of integers
    if not isinstance(seq_pairs, list) or not all(isinstance(pair, list) and len(pair) == 2 and
                                                  all(isinstance(num, int) for num in pair) for pair in seq_pairs):
        raise ValueError(f"Seq_pairs '{seq_pairs}' is not a valid list of [max_seq, extra_seq] pairs.")

    # Check if seeds is an integer
    if not isinstance(seeds, int):
        raise ValueError(f"Seeds '{seeds}' is not a valid integer.")


def run_ensemble_prediction(msa_path, output_path, jobname, max_seq, extra_seq, env, seeds=10, save_all=False):
    complete_output_path = f'{output_path}/{jobname}/predictions/alphafold2/{jobname}'

    if save_all:
        command = (f"colabfold_batch {msa_path} {complete_output_path}_{max_seq}_{extra_seq} "
                   f"--use-dropout --num-seeds {seeds} --max-msa {max_seq}:{extra_seq} --save-all")
    else:
        command = (f"colabfold_batch {msa_path} {complete_output_path}_{max_seq}_{extra_seq} "
                   f"--use-dropout --num-seeds {seeds} --max-msa {max_seq}:{extra_seq}")

    print(
        f'Making prediction for {jobname}'
        f' using max_seq: {max_seq} and extra_seq: {extra_seq} with {seeds} seeds...')

    print(f'\nRead {complete_output_path}_{max_seq}_{extra_seq}/log.txt to monitor output')

    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)

    # Read line by line as they are output
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()
    return rc



def run_multiple_prediction(msa_path, output_path, jobname, seq_pairs, env, seeds=10, save_all=False):
    validate_inputs(msa_path, output_path, jobname, seq_pairs, seeds)
    for max_seq_v, extra_seq_v in seq_pairs:
        run_ensemble_prediction(msa_path,
                                output_path,
                                jobname,
                                max_seq_v,
                                extra_seq_v,
                                env,
                                seeds=seeds, save_all=save_all)


def subset_msa(input_file, output_path, X):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    seqs = []
    seq_name = None
    seq_data = None

    for line in lines:
        if line.startswith('>'):
            if seq_name and seq_data:
                seqs.append((seq_name, seq_data))
            seq_name = line.strip()
            seq_data = ''
        else:
            seq_data += line.strip()

    # Add the last sequence
    if seq_name and seq_data:
        seqs.append((seq_name, seq_data))

    subset_seqs = seqs[:X]

    with open(f'{output_path}/temp_msa.a3m', 'w') as outfile:
        for name, seq in subset_seqs:
            outfile.write(name + '\n')
            outfile.write(seq + '\n')


def load_config(config_file):
    # Default configuration values
    default_config = {
        'msa_path': None,
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

def run_ensemble_prediction(config):
    msa_path = config['msa_path']
    output_path = config['output_path']
    jobname = config['jobname']
    seq_pairs = config['seq_pairs']
    seeds = config['seeds']
    save_all = config['save_all']
    platform = config['platform']
    subset_msa_to = config['subset_msa_to']
    msa_from = config['msa_from']

    if msa_path == None:
        msa_path = f'{output_path}/{jobname}/msas/{msa_from}/{jobname}.a3m'

    print('\n')

    create_directory(f'{output_path}/{jobname}/predictions/alphafold2')

    if subset_msa_to:
        subset_msa(msa_path, output_path, subset_msa_to)

    if platform == 'cpu':
        os.environ['JAX_PLATFORMS'] = platform

    env = os.environ.copy()
    env["PATH"] += os.pathsep + "localcolabfold/colabfold-conda/bin"
    env["PATH"] += os.pathsep + "../localcolabfold/colabfold-conda/bin"
    env["PATH"] += os.pathsep + "/home/gabriel/localcolabfold/colabfold-conda/bin"  ## TODO: remove path

    # Print out the configurations for debugging
    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Job Name: {jobname}")
    print(f"MSA Path: {msa_path}")
    print(f"Output Path: {output_path}")
    print(f"Sequence Pairs: {seq_pairs}")
    print(f"Seeds: {seeds}")
    print(f"Save All: {save_all}")
    if subset_msa_to:
        print(f"Subset MSA to: {subset_msa_to}")
    print(f"MSA From: {msa_from}")
    print("***************************************************************\n")

    # subsets MSA to X sequences
    if subset_msa_to:
        print(f"Subsetting MSA to: {subset_msa_to} sequences\n")
        msa_path = f'{output_path}/{jobname}/predictions/alphafold2/temp_msa.a3m'

    # runs predictions
    run_multiple_prediction(msa_path, output_path, jobname, seq_pairs, env, seeds, save_all)

    # removes temporary MSA if subset
    if subset_msa_to:
        os.remove(f"{output_path}/{jobname}/predictions/alphafold2/temp_msa.a3m")

def main():
    parser = argparse.ArgumentParser(description="Run multiple ensemble predictions.")

    parser.add_argument('--config_file', type=str, help="Path to the configuration file")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--msa_path', type=str, help="Path to the .a3m MSA from jackhmmer")
    parser.add_argument('--output_path', type=str, help="Path to save results to")
    parser.add_argument('--seq_pairs', type=str, help="List of [max_seq, extra_seq] pairs in the format"
                                                      "[[max_seq1, extra_seq1], [max_seq2, extra_seq2], ...]")
    parser.add_argument('--seeds', type=int, nargs='+', help="Number of predictions to run")
    parser.add_argument('--save_all', action='store_true', help="Flag to save all results")
    parser.add_argument('--platform', type=str, choices=['cpu', 'gpu'],
                        help="Platform to run the predictions, GPU recommended")
    parser.add_argument('--subset_msa_to', type=int,
                        help="Subsets the input MSA to X sequences, use when running out of RAM due to deep MSAs")
    parser.add_argument('--msa_from', type=str, choices=['jackhmmer', 'mmseqs2'],
                        help="MSA building tool used to build the input MSA, used to automatically find the MSA if "
                             "--msa_path is not provided")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    config['msa_path'] = args.msa_path if args.msa_path else config.get('msa_path')
    config['output_path'] = args.output_path if args.output_path else config.get('output_path')
    config['jobname'] = args.jobname if args.jobname else config.get('jobname')
    config['seq_pairs'] = args.seq_pairs if args.seq_pairs else config.get('seq_pairs')
    config['seeds'] = args.seeds if args.seeds else config.get('seeds', 10)
    config['save_all'] = args.save_all if 'save_all' in args and args.save_all else config.get('save_all', False)
    config['platform'] = args.platform if args.platform else config.get('platform', 'cpu')
    config['subset_msa_to'] = args.subset_msa_to if args.subset_msa_to else config.get('subset_msa_to')
    config['msa_from'] = args.msa_from if args.msa_from else config.get('msa_from')

    run_ensemble_prediction(config)

if __name__ == "__main__":
    main()
