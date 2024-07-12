import os
import warnings
from decaf_e_dev.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config
from decaf_e_dev.ensemble_analysis.rmsd import rmsd_mode_analysis, build_dataset_rmsd_modes

warnings.filterwarnings("ignore")

def run_rmsd_analysis(config, plot_widget):

    jobname = config.get('jobname')
    output_path = config.get('output_path')
    seq_pairs = config.get('seq_pairs')
    predictions_path = config.get('predictions_path')
    engine = config.get('engine')
    align_range = config.get('align_range')
    analysis_range = config.get('analysis_range')
    analysis_range_name = config.get('analysis_range_name')
    ref1d = config.get('ref1d')
    if ref1d=="":
        ref1d=None
    starting_residue = config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    create_directory(f'{output_path}/{jobname}/analysis/rmsd_1d')
    create_directory(f'{output_path}/{jobname}/analysis/representative_structures/rmsd_1d')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Output Path: {output_path}")
    print(f"max_seq:extra_seq Pairs: {seq_pairs}")
    print(f"Predictions Path: {predictions_path}")
    print(f"Engine: {engine}")
    print(f"Job Name: {jobname}")
    print(f"Align Range: {align_range}")
    print(f"Analysis Range: {analysis_range_name} = {analysis_range}")
    if ref1d:
        print(f"Reference Structure: {ref1d}")
    else:
        print(f"Reference Structure: Top Ranked Prediction")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    print("***************************************************************\n")

    input_dict = {
        'jobname': jobname,
        'output_path': output_path,
        'seq_pairs': seq_pairs,
        'analysis_range': analysis_range,
        'analysis_range_name': analysis_range_name,
        'align_range': align_range
    }

    # Load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # Run 1D RMSD analysis
    rmsd_mode_analysis_dict = rmsd_mode_analysis(pre_analysis_dict, input_dict, ref1d, plot_widget)

    # Build and save results dataset
    build_dataset_rmsd_modes(rmsd_mode_analysis_dict, input_dict)
