import subprocess
import os
import argparse
import json

from fast_ensemble.msa_generation.msa_utils import create_directory


def validate_inputs(config):
    """
    Validate the input configuration parameters.

    Args:
        config (dict): Configuration dictionary containing all necessary parameters.

    Raises:
        ValueError: If any of the configuration parameters are invalid.
    """
    if not os.path.isfile(config['msa_path']):
        raise ValueError(f"MSA path '{config['msa_path']}' is not a valid file.")

    if not os.path.isdir(config['output_path']):
        raise ValueError(f"Output path '{config['output_path']}' is not a valid directory.")

    if not isinstance(config['jobname'], str):
        raise ValueError(f"Jobname '{config['jobname']}' is not a valid string.")

    if not isinstance(config['seq_pairs'], list) or not all(isinstance(pair, list) and len(pair) == 2 and
                                                            all(isinstance(num, int) for num in pair) for pair in config['seq_pairs']):
        raise ValueError(f"Seq_pairs '{config['seq_pairs']}' is not a valid list of [max_seq, extra_seq] pairs.")

    if not isinstance(config['seeds'], int):
        raise ValueError(f"Seeds '{config['seeds']}' is not a valid integer.")


def run_ensemble_prediction_single(msa_path, output_path, jobname, max_seq, extra_seq, env, seeds=10, save_all=False):
    """
    Run a single ensemble prediction with the specified parameters.

    Args:
        msa_path (str): Path to the MSA file.
        output_path (str): Path to save the prediction results.
        jobname (str): Name of the job.
        max_seq (int): Maximum number of sequences to consider in the MSA.
        extra_seq (int): Number of extra sequences to consider in the MSA.
        env (dict): Environment variables for the subprocess.
        seeds (int, optional): Number of predictions to run. Default is 10.
        save_all (bool, optional): Flag to save all results. Default is False.

    Returns:
        int: The return code of the subprocess.
    """
    complete_output_path = f'{output_path}/{jobname}/predictions/alphafold2/{jobname}'

    if save_all:
        command = (f"colabfold_batch {msa_path} {complete_output_path}_{max_seq}_{extra_seq} "
                   f"--use-dropout --num-seeds {seeds} --max-msa {max_seq}:{extra_seq} --save-all")
    else:
        command = (f"colabfold_batch {msa_path} {complete_output_path}_{max_seq}_{extra_seq} "
                   f"--use-dropout --num-seeds {seeds} --max-msa {max_seq}:{extra_seq}")

    print(f'Making prediction for {jobname} using max_seq: {max_seq} and extra_seq: {extra_seq} with {seeds} seeds...')
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


def run_ensemble_prediction(config):
    """
    Run the ensemble predictions based on the provided configuration.

    Args:
        config (dict): Configuration dictionary containing all necessary parameters.
    """
    validate_inputs(config)

    msa_path = config['msa_path']
    output_path = config['output_path']
    jobname = config['jobname']
    seq_pairs = config['seq_pairs']
    seeds = config['seeds']
    save_all = config['save_all']
    platform = config['platform']
    subset_msa_to = config.get('subset_msa_to')
    msa_from = config.get('msa_from', 'jackhmmer')

    if msa_path is None:
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
    env["PATH"] += os.pathsep + "/home/gabriel/localcolabfold/colabfold-conda/bin"  # TODO: remove path

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

    # Subset MSA to X sequences
    if subset_msa_to:
        print(f"Subsetting MSA to: {subset_msa_to} sequences\n")
        msa_path = f"{output_path}/{jobname}/predictions/alphafold2/temp_msa.a3m"

    # Run predictions
    for max_seq_v, extra_seq_v in seq_pairs:
        run_ensemble_prediction_single(msa_path, output_path, jobname, max_seq_v, extra_seq_v, env, seeds, save_all)

    # Remove temporary MSA if subset
    if subset_msa_to:
        os.remove(f"{output_path}/{jobname}/predictions/alphafold2/temp_msa.a3m")


def subset_msa(input_file, output_path, X):
    """
    Subset the MSA file to a specified number of sequences.

    Args:
        input_file (str): Path to the input MSA file.
        output_path (str): Directory path to save the subsetted MSA file.
        X (int): Number of sequences to include in the subset.
    """
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
    """
    Load the configuration from a JSON file.

    Args:
        config_file (str): Path to the configuration file.

    Returns:
        dict: Configuration dictionary with default values if the file is not found or if there's an error in reading the file.
    """
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
    """
    Save the configuration to a JSON file.

    Args:
        config (dict): Configuration dictionary to save.
        file_path (str): Path to save the configuration file.
    """
    with open(file_path, 'w') as file:
        json.dump(config, file, indent=4)


def main():
    """
    Main function to parse command-line arguments and run ensemble predictions.
    """
    parser = argparse.ArgumentParser(description="Run multiple ensemble predictions.")

    parser.add_argument('--config_file', type=str, help="Path to the configuration file")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--msa_path', type=str, help="Path to the .a3m MSA from jackhmmer")
    parser.add_argument('--output_path', type=str, help="Path to save results to")
    parser.add_argument('--seq_pairs', type=str, help="List of [max_seq, extra_seq] pairs in the format [[max_seq1, extra_seq1], [max_seq2, extra_seq2], ...]")
    parser.add_argument('--seeds', type=int, nargs='+', help="Number of predictions to run")
    parser.add_argument('--save_all', action='store_true', help="Flag to save all results")
    parser.add_argument('--platform', type=str, choices=['cpu', 'gpu'], help="Platform to run the predictions, GPU recommended")
    parser.add_argument('--subset_msa_to', type=int, help="Subsets the input MSA to X sequences, use when running out of RAM due to deep MSAs")
    parser.add_argument('--msa_from', type=str, choices=['jackhmmer', 'mmseqs2'], help="MSA building tool used to build the input MSA, used to automatically find the MSA if --msa_path is not provided")

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