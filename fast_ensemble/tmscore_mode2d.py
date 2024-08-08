import os
import warnings
from fast_ensemble.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config, auto_select_2d_references
from fast_ensemble.ensemble_analysis.twotmscore import TwoTMScore

warnings.filterwarnings("ignore")

def run_2d_tmscore_analysis(config, widget):

    output_path = config.get('output_path')
    predictions_path = config.get('predictions_path')
    mode_results = config.get('mode_results')
    seq_pairs = config.get('seq_pairs')
    jobname = config.get('jobname')
    ref2d1 = config.get('ref2d1')
    ref2d2 = config.get('ref2d2')
    starting_residue = config.get('starting_residue')
    slice_predictions = config.get('slice_predictions')
    engine = config.get('engine')
    n_stdevs = config.get('n_stdevs')
    n_clusters = config.get('n_clusters')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    if not mode_results:
        mode_results = f'{output_path}/{jobname}/analysis/tmscore_1d/{jobname}_tmscore_1d_analysis_results.csv'

    if not ref2d1 and not ref2d2:
        ref2d1, ref2d2 = auto_select_2d_references(mode_results, 'tmscore')

    create_directory(f'{output_path}/{jobname}/analysis/tmscore_2d')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Predictions Path: {predictions_path}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")
    print(f"Engine: {engine}")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    if slice_predictions:
        print(f"Setting Analysis Range to: {slice_predictions}")
    print(f"Reference 1: {ref2d1}")
    print(f"Reference 2: {ref2d2}")
    print(f"Number of Standard Devs. to Consider Point Closeness: {n_stdevs}")
    if n_clusters:
        print(f"Number of Clusters: {n_clusters}")
    else:
        print(f"Number of Clusters: Number of Detected 1D TM-Score Modes + 1")
    print("***************************************************************\n")

    input_dict = {'jobname': jobname,
                  'output_path': output_path,
                  'seq_pairs': seq_pairs,
                  'predictions_path': predictions_path}

    # loads predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # runs 2D TM-Score analysis
    twod = TwoTMScore(pre_analysis_dict, input_dict, widget, ref2d1, ref2d2, slice_predictions)

    # builds results dataset and saves to disk
    twod.get_2d_tmscore(mode_results, n_stdevs, n_clusters, output_path)
