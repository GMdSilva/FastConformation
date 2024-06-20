import pandas as pd

from tqdm import tqdm

from ensemble_analysis.analysis_utils import save_traj, reorder_frames_by


TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'


def save_trajs(prediction_dicts, input_dict, reorder, traj_format):
    jobname = input_dict['jobname']
    analysis_range = input_dict['analysis_range']
    analysis_range_name = input_dict['analysis_range_name']
    output_path = input_dict['output_path']

    full_output_path = f"{output_path}/{jobname}/trajs/"
    print('')
    with tqdm(total=len(prediction_dicts), bar_format=TQDM_BAR_FORMAT) as pbar:
        for prediction in prediction_dicts:
            pbar.set_description(f'Saving trajectory for {prediction}')
            frames = prediction_dicts[prediction]['mda_universe']
            max_seq = prediction_dicts[prediction]['max_seq']
            extra_seq = prediction_dicts[prediction]['extra_seq']

            if reorder:
                order = pd.read_csv(f'{output_path}/'
                                    f'{jobname}/'
                                    f'analysis/'
                                    f'{reorder}/'
                                    f'{jobname}_'
                                    f'{max_seq}_'
                                    f'{extra_seq}_'
                                    f'{reorder}_'
                                    f'df.csv')

                if reorder == 'pca':
                    frames = reorder_frames_by(frames, order['PC1'])
                elif reorder == 'tmscore':
                    frames = reorder_frames_by(frames, order['tmscore'])
                else:
                    frames = reorder_frames_by(frames, order[analysis_range])

            save_traj(frames, full_output_path, jobname, max_seq, extra_seq, traj_format=traj_format, ordered=reorder)
            pbar.update(n=1)

    print(f'\nSaved {jobname} trajectories to {full_output_path}\n')

    return prediction_dicts
