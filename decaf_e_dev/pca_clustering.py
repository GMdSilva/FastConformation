import os
import warnings
from decaf_e_dev.ensemble_analysis.analysis_utils import create_directory, load_predictions, load_config
from decaf_e_dev.ensemble_analysis.pca import pca_from_ensemble

warnings.filterwarnings("ignore")

def run_pca_analysis(config):

    predictions_path = config.get('predictions_path')
    output_path = config.get('output_path')
    seq_pairs = config.get('seq_pairs')
    jobname = config.get('jobname')
    align_range = config.get('align_range')
    analysis_range = config.get('analysis_range')
    analysis_range_name = config.get('analysis_range_name')
    engine = config.get('engine')
    n_pca_clusters = config.get('n_pca_clusters')
    starting_residue = config.get('starting_residue')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f"Output path {output_path} is not a directory")

    if not predictions_path:
        predictions_path = f'{output_path}/{jobname}/predictions/{engine}'

    create_directory(f'{output_path}/{jobname}/analysis/pca')

    print("\nConfigurations:")
    print("***************************************************************")
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

