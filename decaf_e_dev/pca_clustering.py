import os
import warnings
from decaf_e_dev.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config
from decaf_e_dev.ensemble_analysis.pca import pca_from_ensemble

warnings.filterwarnings("ignore")

def run_pca_analysis(config_file=None, predictions_path=None, output_path=None, seq_pairs=None,
                     jobname=None, align_range=None, analysis_range=None, analysis_range_name=None,
                     engine=None, n_pca_clusters=None, starting_residue=None):

    # Load configuration from file if provided
    config_file = config_file if config_file else 'config.json'
    config = load_config(config_file)

    # Override config with function arguments if provided
    predictions_path = predictions_path if predictions_path else config.get('predictions_path')
    output_path = output_path if output_path else config.get('output_path')
    seq_pairs = seq_pairs if seq_pairs else config.get('seq_pairs')
    jobname = jobname if jobname else config.get('jobname')
    align_range = align_range if align_range else config.get('align_range')
    analysis_range = analysis_range if analysis_range else config.get('analysis_range')
    analysis_range_name = analysis_range_name if analysis_range_name else config.get('analysis_range_name')
    engine = engine if engine else config.get('engine')
    n_pca_clusters = n_pca_clusters if n_pca_clusters else config.get('n_pca_clusters')
    starting_residue = starting_residue if starting_residue else config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    create_directory(f'{output_path}/{jobname}/analysis/pca')

    print("\nConfigurations:")
    print("***************************************************************")
    print(f"Used Config File? {config_file is not None}")
    print(f"Predictions Path: {predictions_path}")
    print(f"Output Path: {output_path}")
    print(f"Job Name: {jobname}")
    print(f"Analysis Range: {analysis_range_name} = {analysis_range}")
    print(f"Align Range: {align_range}")
    print(f"Number of Clusters: {n_pca_clusters}")
    if starting_residue:
        print(f"Starting Residue: {starting_residue}")
    print("***************************************************************\n")

    # load predictions to RAM
    pre_analysis_dict = load_predictions(predictions_path, seq_pairs, jobname, starting_residue)

    # run PCA analysis
    pca_from_ensemble(jobname, pre_analysis_dict, output_path, align_range, analysis_range, int(n_pca_clusters))

