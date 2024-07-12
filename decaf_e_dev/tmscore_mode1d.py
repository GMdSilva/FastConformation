import os
import warnings
from decaf_e_dev.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config
from decaf_e_dev.ensemble_analysis.tmscore import build_dataset_tmscore_modes, tmscore_mode_analysis

warnings.filterwarnings("ignore")

def run_tmscore_analysis(config, plot_widget):

    output_path = config.get('output_path')
    predictions_path = config.get('predictions_path')
    seq_pairs = config.get('seq_pairs')
    jobname = config.get('jobname')
    slice_predictions = config.get('slice_predictions')
    ref1 = config.get('ref1')
    engine = config.get('engine')
    starting_residue = config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    create_directory(f'{output_path}/{jobname}/analysis/tmscore_1d')
    create_directory(f'{output_path}/{jobname}/analysis/representative_structures/tmscore_1d')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Predictions Path: {predictions_path}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")
    print(f"Engine: {engine}")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    if ref1:
        print(f"Reference: {ref1}")
    if slice_predictions:
        print(f"Setting Analysis Range to: {slice_predictions}")
    print("***************************************************************\n")

    input_dict = {
        'jobname': jobname,
        'predictions_path': predictions_path,
        'output_path': output_path,
        'seq_pairs': seq_pairs
    }

    # Load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)
    print("Load predictions to RAM")
    # Run 1D TM-Score analysis
    tmscore_mode_analysis_dict = tmscore_mode_analysis(pre_analysis_dict, input_dict, ref1, slice_predictions, plot_widget)
    print("Run 1D TM-Score analysis")
    # Build results dataset and save to disk
    build_dataset_tmscore_modes(tmscore_mode_analysis_dict, input_dict)
