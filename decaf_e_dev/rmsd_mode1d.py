import os
import warnings
from decaf_e_dev.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config
from decaf_e_dev.ensemble_analysis.rmsd import rmsd_mode_analysis, build_dataset_rmsd_modes

warnings.filterwarnings("ignore")

def run_rmsd_analysis(config_file=None, jobname=None, output_path=None, seq_pairs=None,
                      predictions_path=None, engine=None, align_range=None,
                      analysis_range=None, analysis_range_name=None,
                      ref1d=None, starting_residue=None):

    # Load configuration from file if provided
    config_file = config_file if config_file else 'config.json'
    config = load_config(config_file)

    # Override config with function arguments if provided
    jobname = jobname if jobname else config.get('jobname')
    output_path = output_path if output_path else config.get('output_path')
    seq_pairs = seq_pairs if seq_pairs else config.get('seq_pairs')
    predictions_path = predictions_path if predictions_path else config.get('predictions_path')
    engine = engine if engine else config.get('engine')
    align_range = align_range if align_range else config.get('align_range')
    analysis_range = analysis_range if analysis_range else config.get('analysis_range')
    analysis_range_name = analysis_range_name if analysis_range_name else config.get('analysis_range_name')
    ref1d = ref1d if ref1d else config.get('ref1d')
    starting_residue = starting_residue if starting_residue else config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    create_directory(f'{output_path}/{jobname}/analysis/rmsd_1d')
    create_directory(f'{output_path}/{jobname}/analysis/representative_structures/rmsd_1d')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Used Config File? {config_file is not None}")
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
    rmsd_mode_analysis_dict = rmsd_mode_analysis(pre_analysis_dict, input_dict, ref1d)

    # Build and save results dataset
    build_dataset_rmsd_modes(rmsd_mode_analysis_dict, input_dict)
