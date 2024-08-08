import os
import warnings
from fast_ensemble.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config, auto_select_2d_references
from fast_ensemble.ensemble_analysis.twodrmsd import TwodRMSD

warnings.filterwarnings("ignore")

def run_2d_rmsd_analysis(config, widget):

    jobname = config.get('jobname')
    output_path = config.get('output_path')
    mode_results = config.get('mode_results')
    seq_pairs = config.get('seq_pairs')
    predictions_path = config.get('predictions_path')
    engine = config.get('engine')
    align_range = config.get('align_range')
    analysis_range = config.get('analysis_range')
    analysis_range_name = config.get('analysis_range_name')
    ref2d1 = config.get('ref2d1')
    ref2d2 = config.get('ref2d2')
    n_stdevs = config.get('n_stdevs')
    n_clusters = config.get('n_clusters')
    starting_residue = config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    if not mode_results:
        mode_results = f'{output_path}/{jobname}/analysis/rmsd_1d/{jobname}_rmsd_1d_analysis_results.csv'

    if not ref2d1 and not ref2d2:
        ref2d1, ref2d2 = auto_select_2d_references(mode_results, 'RMSD')

    create_directory(f'{output_path}/{jobname}/analysis/rmsd_2d')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Job Name: {jobname}")
    print(f"Output Path: {output_path}")
    print(f"max_seq:extra_seq Pairs: {seq_pairs}")
    print(f"Predictions Path: {predictions_path}")
    print(f"Engine: {engine}")
    print(f"Align Range: {align_range}")
    print(f"Analysis Range: {analysis_range_name} = {analysis_range}")
    print(f"Reference 1: {ref2d1}")
    print(f"Reference 2: {ref2d2}")
    print(f"Number of Standard Devs. to Consider Points Close to fit Line: {n_stdevs}")
    if n_clusters:
        print(f"Number of Clusters: {n_clusters}")
    else:
        print(f"Number of Clusters: Number of Detected 1D RMSD Modes + 1")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    print("***************************************************************\n")

    input_dict = {
        'jobname': jobname,
        'output_path': output_path,
        'seq_pairs': seq_pairs,
        'analysis_range': analysis_range,
        'analysis_range_name': analysis_range_name,
        'align_range': align_range,
    }

    # Load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # Run 2D RMSD analysis
    twod = TwodRMSD(pre_analysis_dict, input_dict, widget, ref2d1, ref2d2)

    # Build and save results dataset
    twod.get_2d_rmsd(mode_results, n_stdevs, n_clusters, output_path)
