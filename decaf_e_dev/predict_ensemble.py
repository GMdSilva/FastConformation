import subprocess
import os
import argparse
import json


def validate_inputs(msa_path, output_path, jobname, seq_pairs, seeds):
    # Check if msa_path is a valid file
    if not os.path.isfile(msa_path):
        raise ValueError(f"MSA path '{msa_path}' is not a valid file.")

    if not os.path.isdir(output_path):
        raise ValueError(f"Output path '{msa_path}' is not a valid directory.")

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
    complete_output_path = f'{output_path}/{jobname}'

    if save_all:
        command = (f"colabfold_batch {msa_path} {complete_output_path}_preds_{max_seq}_{extra_seq}/"
                   f"--use-dropout --num-seeds {seeds} --max-msa {max_seq}:{extra_seq} --save-all")
    else:
        command = (f"colabfold_batch {msa_path} {complete_output_path}_preds_{max_seq}_{extra_seq}/"
                   f"--use-dropout --num-seeds {seeds} --max-msa {max_seq}:{extra_seq}")

    print(
        f'Making prediction for {jobname}'
        f' using max_seq: {max_seq} and extra_seq: {extra_seq} with {seeds} seeds...')

    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                               env=env)

    # Read line by line as they are output
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()


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


def main():
    parser = argparse.ArgumentParser(description="Run multiple ensemble predictions.")

    parser.add_argument('--config_file', type=str, default='config.json',
                        help="Path to the configuration file (default: config.json)")
    parser.add_argument('--msa_path', type=str, help="Path to the .a3m MSA from jackhmmer")
    parser.add_argument('--output_path', type=str, help="Path to save results to"
                                                        "(default: current directory)")
    parser.add_argument('--jobname', type=str, help="The job name")
    parser.add_argument('--seq_pairs', type=str, help="List of [max_seq, extra_seq] pairs in the format"
                                                      "[[max_seq1, extra_seq1], [max_seq2, extra_seq2], ...]")
    parser.add_argument('--seeds', type=int, nargs='+', help="Number of predictions to run")
    parser.add_argument('--save_all', action='store_true', help="Flag to save all results")
    parser.add_argument('--platform', type=str, choices=['cpu', 'gpu'],
                        help="Platform to run the predictions (default: cpu)")

    args = parser.parse_args()

    # Load configuration from file if provided
    config_file = args.config_file if args.config_file else 'config.json'
    config = load_config(config_file)

    # Override config with command line arguments if provided
    msa_path = args.msa_path if args.msa_path else config.get('msa_path')
    output_path = args.output_path if args.output_path else config.get('output_path')
    jobname = args.jobname if args.jobname else config.get('jobname')
    seq_pairs = args.seq_pairs if args.seq_pairs else config.get('seq_pairs')
    seeds = args.seeds if args.seeds else config.get('seeds', 10)
    save_all = args.save_all if 'save_all' in args and args.save_all else config.get('save_all', False)
    platform = args.platform if args.platform else config.get('platform', 'cpu')

    if platform == 'cpu':
        os.environ['JAX_PLATFORMS'] = platform

    env = os.environ.copy()
    env["PATH"] += os.pathsep + "localcolabfold/colabfold-conda/bin"
    env["PATH"] += os.pathsep + "/home/gabriel/localcolabfold/colabfold-conda/bin"  ## TODO: remove path

    # Print out the configurations for debugging
    print("Configurations:")
    print(f"Used Config File? {args.config_file is not None}")
    print(f"MSA Path: {msa_path}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")
    print(f"Sequence Pairs: {seq_pairs}")
    print(f"Seeds: {seeds}")
    print(f"Save All: {save_all}")
    print(f"Platform: {platform}")

    run_multiple_prediction(msa_path, output_path, jobname, seq_pairs, env, seeds, save_all)

    save_config(config, config_file)


if __name__ == "__main__":
    main()
