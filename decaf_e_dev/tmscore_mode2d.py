import os
import warnings
from decaf_e_dev.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config, auto_select_2d_references
from decaf_e_dev.ensemble_analysis.twotmscore import TwoTMScore

warnings.filterwarnings("ignore")

def run_2d_tmscore_analysis(config_file=None, output_path=None, predictions_path=None, mode_results=None, jobname=None, seq_pairs=None,
                            starting_residue=None, slice_predictions=None, ref2d1=None, ref2d2=None, engine=None, n_stdevs=None, n_clusters=None):

    # Load configuration from file if provided
    config_file = config_file if config_file else 'config.json'
    config = load_config(config_file)

    # Override config with function arguments if provided
    output_path = output_path if output_path else config.get('output_path')
    predictions_path = predictions_path if predictions_path else config.get('predictions_path')
    mode_results = mode_results if mode_results else config.get('mode_results')
    seq_pairs = seq_pairs if seq_pairs else config.get('seq_pairs')
    jobname = jobname if jobname else config.get('jobname')
    ref2d1 = ref2d1 if ref2d1 else config.get('ref2d1')
    ref2d2 = ref2d2 if ref2d2 else config.get('ref2d2')
    starting_residue = starting_residue if starting_residue else config.get('starting_residue')
    slice_predictions = slice_predictions if slice_predictions else config.get('slice_predictions')
    engine = engine if engine else config.get('engine')
    n_stdevs = n_stdevs if n_stdevs else config.get('n_stdevs')
    n_clusters = n_clusters if n_clusters else config.get('n_clusters')

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
    print(f"Used Config File? {config_file is not None}")
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
    twod = TwoTMScore(pre_analysis_dict, input_dict, ref2d1, ref2d2, slice_predictions)

    # builds results dataset and saves to disk
    twod.get_2d_tmscore(mode_results, n_stdevs, n_clusters)
