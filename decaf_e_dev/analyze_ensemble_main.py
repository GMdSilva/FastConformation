from ensemble_analysis.analysis_utils import load_config
from ensemble_analysis.analysis_utils import create_directory
from ensemble_analysis.analysis_utils import load_predictions, load_predictions_json
from ensemble_analysis.rmsd import rmsd_mode_analysis, build_dataset_rmsd_modes
from ensemble_analysis.twodrmsd import TwodRMSD
#from scripts.rmsf_analysis import calculate_rmsf_multiple, calculate_rmsf_and_call_peaks, build_dataset_rmsf_peaks
#from scripts.pca_analysis_script import pca_from_ensemble_movie

config = load_config('example_configs/config.json')

# Override config with command line arguments if provided
output_path = config.get('output_path')
seq_pairs = config.get('seq_pairs')
jobname = config.get('jobname')
align_range = config.get('align_range')
analysis_range = config.get('analysis_range')
analysis_range_name = config.get('analysis_range_name')
#
# create_directory(f'{output_path}/representative_structures')
# create_directory(f'{output_path}/plots')
# create_directory(f'{output_path}/datasets')

# rmsd_mode_analysis_dict = rmsd_mode_analysis(pre_analysis_dict, align_range, rmsd_ranges, output_path)
# build_dataset_rmsd_modes(rmsd_mode_analysis_dict, rmsd_ranges, output_path)

input_dict = {'jobname': jobname,
              'output_path': output_path,
              'seq_pairs': seq_pairs,
              'rmsd_range': analysis_range,
              'rmsd_range_name': analysis_range_name,
              'align_range': align_range,}
pre_analysis_dict = load_predictions(output_path, seq_pairs, jobname)
#rmsd_mode_analysis_dict = rmsd_mode_analysis(pre_analysis_dict, input_dict)
# aa = build_dataset_rmsd_modes(rmsd_mode_analysis_dict, input_dict)


twod = TwodRMSD(pre_analysis_dict, input_dict)
twod.get_2d_rmsd('sample_predictions/datasets/abl_wt_rmsd_mode_analysis_results.csv')



#calculate_rmsf_multiple(jobname, pre_analysis_dict, align_range, output_path)
#rmsf_peak_calling_dict = calculate_rmsf_and_call_peaks(jobname, pre_analysis_dict, align_range, output_path)
#build_dataset_rmsf_peaks(rmsf_peak_calling_dict, output_path)

#rmsd_ranges = [f'backbone and resid {peak["starting_residue"]}-{peak["ending_residue"]}' for peak in rmsf_peak_calling_dict]
#
# rmsd_mode_analysis_dict = rmsd_mode_analysis(pre_analysis_dict, align_range, rmsd_ranges, output_path)
# build_dataset_rmsd_modes(rmsd_mode_analysis_dict, rmsd_ranges, output_path)
#
# pca_from_ensemble_movie(jobname, pre_analysis_dict, output_path, align_range, analysis_ranges)