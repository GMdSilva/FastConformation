import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from scipy.optimize import curve_fit

from tqdm import tqdm

from decaf_e_dev.ensemble_analysis.analysis_utils import parabola
from decaf_e_dev.ensemble_analysis.rmsd import calculate_rmsd


TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'


class TwodRMSD:
    def __init__(self, prediction_dicts, input_dict, ref_gr=None, ref_alt=None):
        self.prediction_dicts = prediction_dicts
        self.input_dict = input_dict
        self.ref_gr = ref_gr
        self.ref_alt = ref_alt
        self.filtering_dict = {}
        self.clustering_dict = {}

    def calculate_2d_rmsd(self, trial):
        universe = self.prediction_dicts[trial]['mda_universe']
        rmsd_gr = calculate_rmsd(universe,
                                 ref=self.ref_gr,
                                 align_range='backbone',
                                 analysis_range=self.input_dict['analysis_range'])

        rmsd_alt = calculate_rmsd(universe,
                                  ref=self.ref_alt,
                                  align_range='backbone',
                                  analysis_range=self.input_dict['analysis_range'])
        rmsd_2d_data = np.array([rmsd_gr[self.input_dict['analysis_range']],
                                 rmsd_alt[self.input_dict['analysis_range']]]).T

        return rmsd_2d_data

    def fit_and_filter_data(self, rmsd_2d_data, n_stdevs):
        popt, _ = curve_fit(parabola, rmsd_2d_data[:, 0], rmsd_2d_data[:, 1])
        fit_x = np.linspace(min(rmsd_2d_data[:, 0]), max(rmsd_2d_data[:, 0]), 100)
        fit_y = parabola(fit_x, *popt)

        fitted_curve_values = parabola(rmsd_2d_data[:, 0], *popt)
        distances = np.abs(rmsd_2d_data[:, 1] - fitted_curve_values)

        threshold = np.mean(distances) + int(n_stdevs) * np.std(distances)

        close_points = distances < threshold
        x_close = rmsd_2d_data[close_points, 0]
        y_close = rmsd_2d_data[close_points, 1]

        bin_edges = np.linspace(min(x_close), max(x_close), 101)
        bins = np.digitize(x_close, bin_edges)

        unique_bins = np.unique(bins)
        ratio = len(unique_bins) / 100

        self.filtering_dict = {
            'fit_x': fit_x,
            'fit_y': fit_y,
            'close_points': close_points,
            'bins': bins,
            'unique_bins': unique_bins,
            'ratio': ratio,
            'bin_edges': bin_edges,
            'x_close': x_close,
            'y_close': y_close,
        }

    def plot_filtering_data(self, rmsd_2d_data):
        plt.figure(figsize=(5, 4))
        plt.scatter(rmsd_2d_data[:, 0], rmsd_2d_data[:, 1], s=10)

        plt.plot(self.filtering_dict['fit_x'],
                 self.filtering_dict['fit_y'],
                 label='Fitted Curve', color='red')

        plt.scatter(self.filtering_dict['x_close'],
                    self.filtering_dict['y_close'],
                    label='Close Points',
                    color='green',
                    s=20)

        title = (f"{self.input_dict['jobname']} "
                 f"{self.input_dict['max_seq']} "
                 f"{self.input_dict['extra_seq']}")

        plt.title(title, fontsize=15)
        plt.xlabel('RMSD vs. Ref1 (A)', fontsize=14)
        plt.ylabel('RMSD vs. Ref2 (A)', fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.legend(loc='best')
        plt.tight_layout()

        plot_path = (f"{self.input_dict['output_path']}/"
                     f"{self.input_dict['jobname']}/"
                     f"analysis/"
                     f"rmsd_2d/"
                     f"{self.input_dict['jobname']}_"
                     f"{self.input_dict['max_seq']}_"
                     f"{self.input_dict['extra_seq']}_"
                     f"rmsd_2d_fit.png")

        plt.savefig(plot_path, dpi=300)
        plt.close()

    def cluster_2d_data(self, rmsd_2d_data, n_clusters):
        kmeans = KMeans(n_clusters)
        close_points_2d = np.array([self.filtering_dict['x_close'],
                                    self.filtering_dict['y_close']]).T

        kmeans.fit(close_points_2d)
        labels = kmeans.labels_
        centroids = kmeans.cluster_centers_
        sorted_indices = np.argsort(centroids[:, 0])
        correct_labels = np.zeros_like(labels)

        for correct_label, original_index in enumerate(sorted_indices):
            correct_labels[labels == original_index] = correct_label

        unique_labels = set(correct_labels)
        cluster_counts = {i: 0 for i in unique_labels}
        total_samples = close_points_2d.shape[0]
        outliers = rmsd_2d_data.shape[0] - total_samples
        outliers = (outliers / total_samples) * 100

        for k in unique_labels:
            class_member_mask = (correct_labels == k)
            xy = close_points_2d[class_member_mask]
            cluster_counts[k] = round((len(xy) / total_samples) * 100, 1)
        self.clustering_dict = {
            'labels': labels,
            'correct_labels': correct_labels,
            'centroids': centroids,
            'cluster_counts': cluster_counts,
            'outliers': outliers,
            'unique_labels': unique_labels,
            'close_points_2d': close_points_2d
        }

    def plot_and_save_2d_data(self):
        unique_labels = self.clustering_dict['unique_labels']
        correct_labels = self.clustering_dict['correct_labels']
        cluster_counts = self.clustering_dict['cluster_counts']
        centroids = self.clustering_dict['centroids']
        outliers = self.clustering_dict['outliers']

        plt.figure(figsize=(5, 4))
        colors = ['blue', 'green', 'magenta', 'orange', 'grey', 'brown', 'cyan', 'purple']

        for i in unique_labels:
            cluster_points = self.clustering_dict['close_points_2d'][correct_labels == i]
            plt.scatter(cluster_points[:, 0], cluster_points[:, 1], c=colors[i],
                        label=f'Cluster {i} pop: {cluster_counts[i]}', alpha=0.6)
        plt.scatter(centroids[:, 0], centroids[:, 1], s=100, c='black', marker='X', label='Centroids')
        plt.plot(self.filtering_dict['fit_x'], self.filtering_dict['fit_y'], label='Fitted Curve', color='red')

        title = (f"{self.input_dict['jobname']} "
                 f"{self.input_dict['max_seq']} "
                 f"{self.input_dict['extra_seq']} "
                 f"Score: {self.filtering_dict['ratio']:.2f}")

        plt.title(title, fontsize=16)
        plt.legend()
        plt.xlabel('RMSD vs. Ref1 (A)', fontsize=14)
        plt.ylabel('RMSD vs. Ref2 (A)', fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.tight_layout()

        plot_path = (f"{self.input_dict['output_path']}/"
                     f"{self.input_dict['jobname']}/"
                     f"analysis/"
                     f"rmsd_2d/"
                     f"{self.input_dict['jobname']}_"
                     f"{self.input_dict['max_seq']}_"
                     f"{self.input_dict['extra_seq']}_"
                     f"rmsd_2d_clustered.png")

        plt.savefig(plot_path, dpi=300)
        plt.close()

        records = []

        for i in unique_labels:
            records.append({
                'trial': self.input_dict['trial'],
                'analysis_range': self.input_dict['analysis_range'],
                'score': self.filtering_dict['ratio'],
                'cluster_label': i,
                'cluster_pop': cluster_counts[i],
                'centroid_values': centroids[i],
                '%_outliers': outliers
            })

        df = pd.DataFrame(records)
        return df

    def get_2d_rmsd(self, rmsd_mode_df_path, n_stdevs, n_clusters):
        df = pd.read_csv(rmsd_mode_df_path)
        unique_trials = df['trial'].unique()

        df_all_trials = pd.DataFrame()
        with tqdm(total=len(self.prediction_dicts), bar_format=TQDM_BAR_FORMAT) as pbar:
            print("jjjjj")
            for trial in unique_trials:
                if not n_clusters:
                    unique_df = df[df['trial'] == trial]
                    n_clusters_trial = len(unique_df['mode_label']) + 1
                pbar.set_description(f'Running 2D RMSD analysis for {trial}')
                self.input_dict['trial'] = trial
                self.input_dict['max_seq'] = self.prediction_dicts[trial]['max_seq']
                self.input_dict['extra_seq'] = self.prediction_dicts[trial]['extra_seq']

                rmsd_2d_data = self.calculate_2d_rmsd(trial)
                if len(rmsd_2d_data) > 0:
                    print("jjj")
                    self.fit_and_filter_data(rmsd_2d_data, n_stdevs)
                    print("jhhhjj")
                    self.plot_filtering_data(rmsd_2d_data)
                    print("jhhhjjkkk")
                    self.cluster_2d_data(rmsd_2d_data, n_clusters_trial)
                    df_to_save = self.plot_and_save_2d_data()
                    print("rrr")
                    df_all_trials = pd.concat([df_all_trials, df_to_save], ignore_index=True)
                    print("rrrrrrerrrrr")
                pbar.update(n=1)

        csv_path = (f"{self.input_dict['output_path']}/"
                    f"{self.input_dict['jobname']}/"
                    f"analysis/"
                    f"rmsd_2d/"
                    f"{self.input_dict['jobname']}_"
                    f"clustering_analysis.csv")

        print(f"\nSaving {self.input_dict['jobname']} 2D RMSD analysis results to {csv_path}\n")
        df_all_trials.to_csv(csv_path, index=False)
