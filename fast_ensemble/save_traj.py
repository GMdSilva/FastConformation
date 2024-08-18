import os
import warnings
from fast_ensemble.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config
from fast_ensemble.ensemble_analysis.traj import save_trajs

warnings.filterwarnings("ignore")

def run_trajectory_saving(config):
    """
    Run the process of saving ensemble predictions as trajectories.

    Args:
        config (dict): Configuration dictionary containing the following keys:
            - output_path (str): Path to save the output trajectories.
            - predictions_path (str, optional): Path to the predictions directory. If not provided, it defaults to the
              predictions directory within the output path.
            - seq_pairs (list): List of sequence pairs for the analysis.
            - jobname (str): Name of the job.
            - analysis_range (tuple): The range of residues or frames to analyze.
            - analysis_range_name (str): The name associated with the analysis range.
            - reorder (str, optional): Reordering method to use (e.g., RMSD, TMSCORE, PCA).
            - traj_format (str): Format to save the trajectory (e.g., PDB, DCD).
            - engine (str): The molecular dynamics engine used (e.g., AlphaFold, Rosetta).
            - starting_residue (int, optional): The starting residue number for the trajectory.

    Raises:
        NotADirectoryError: If the specified output path is not a directory.
    """
    output_path = config.get('output_path')
    predictions_path = config.get('predictions_path')
    seq_pairs = config.get('seq_pairs')
    jobname = config.get('jobname')
    analysis_range = config.get('analysis_range')
    analysis_range_name = config.get('analysis_range_name')
    reorder = config.get('reorder')
    traj_format = config.get('traj_format')
    engine = config.get('engine')
    starting_residue = config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    create_directory(f'{output_path}/{jobname}/trajs')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Predictions Path: {predictions_path}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")
    print(f"Analysis Range: {analysis_range_name} = {analysis_range}")
    print(f"Reorder by: {reorder}")
    print(f"Traj format: {traj_format}")
    print(f"Engine: {engine}")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    print("***************************************************************\n")

    input_dict = {'jobname': jobname,
                  'output_path': output_path,
                  'seq_pairs': seq_pairs,
                  'analysis_range': analysis_range,
                  'analysis_range_name': analysis_range_name}

    # load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # save predictions as a trajectory, optionally reordered by X analysis (RMSD 1D, TMSCORE 1D, PCA, etc.)
    save_trajs(pre_analysis_dict, input_dict, reorder, traj_format)
