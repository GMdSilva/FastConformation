import re
import os
import shutil
import subprocess
from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from scipy.optimize import curve_fit
import MDAnalysis as mda
from tqdm import tqdm
from fast_ensemble.ensemble_analysis.analysis_utils import parabola
from fast_ensemble.ensemble_analysis.analysis_utils import create_directory


TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'

class TwoTMScore:
    """
    A class to perform 2D TM-Score analysis on molecular dynamics simulations.

    Attributes:
    ----------
    prediction_dicts : dict
        A dictionary containing prediction data with associated MDAnalysis Universes.
    input_dict : dict
        A dictionary containing job-related metadata (jobname, analysis range, etc.).
    slice_predictions : str or None
        A selection string for slicing the predictions (default is None).
    ref_gr : str or None
        The reference structure file path for the first TM-Score calculation (default is None).
    ref_alt : str or None
        The reference structure file path for the second TM-Score calculation (default is None).
    filtering_dict : dict
        A dictionary to store data related to the filtering of TM-Score values.
    clustering_dict : dict
        A dictionary to store data related to the clustering of 2D TM-Score values.
    widget : object
        A widget object to handle the plotting of the analysis results.

    Methods:
    -------
    slice_models(universe, selection, temp_traj_path):
        Static method to slice models according to a given selection and save the results.
    tmscore_wrapper(mobile, reference):
        Static method to run the TM-Score command and return the TM-Score value.
    run_tmscore(folder_path, custom_ref):
        Run TM-Score on all PDB files in a folder against a custom reference structure.
    calculate_2d_tmscore(trial):
        Calculate 2D TM-Score for a given trial.
    fit_and_filter_data(tmscore_2d_data, n_stdevs):
        Fit a parabola to the 2D TM-Score data and filter points based on the standard deviation threshold.
    plot_filtering_data(tmscore_2d_data):
        Generate and save a plot of the filtered 2D TM-Score data with the fitted curve.
    cluster_2d_data(tmscore_2d_data, n_clusters):
        Perform clustering on the filtered 2D TM-Score data and store clustering results.
    plot_and_save_2d_data(output_path):
        Plot the clustered 2D TM-Score data, save the plot, and return a DataFrame with the clustering information.
    get_2d_tmscore(tmscore_mode_df_path, n_stdevs, n_clusters, output_path):
        Execute the full 2D TM-Score analysis for all trials, including fitting, filtering, clustering, and saving results.

    """

    def __init__(self, prediction_dicts, input_dict, widget, ref_gr=None, ref_alt=None, slice_predictions=None):
        """
        Initialize the TwoTMScore class.

        Parameters:
        ----------
        prediction_dicts : dict
            A dictionary containing prediction data with associated MDAnalysis Universes.
        input_dict : dict
            A dictionary containing job-related metadata (jobname, analysis range, etc.).
        widget : object
            A widget object to handle the plotting of the analysis results.
        ref_gr : str or None, optional
            The reference structure file path for the first TM-Score calculation (default is None).
        ref_alt : str or None, optional
            The reference structure file path for the second TM-Score calculation (default is None).
        slice_predictions : str or None, optional
            A selection string for slicing the predictions (default is None).
        """
        self.prediction_dicts = prediction_dicts
        self.input_dict = input_dict
        self.slice_predictions = slice_predictions
        self.ref_gr = ref_gr
        self.ref_alt = ref_alt
        self.filtering_dict = {}
        self.clustering_dict = {}
        self.widget = widget

    @staticmethod
    def slice_models(universe, selection, temp_traj_path):
        """
        Static method to slice models according to a given selection and save the results.

        Parameters:
        ----------
        universe : MDAnalysis.Universe
            The Universe containing the trajectory data.
        selection : str
            The selection string for slicing the models.
        temp_traj_path : str
            The path where the sliced models will be saved.

        Returns:
        -------
        None
        """
        for idx, ts in enumerate(universe.trajectory):
            ag = universe.select_atoms(selection)
            ag.write(f"{temp_traj_path}/tmp_{idx:0>5}.pdb")

    @staticmethod
    def tmscore_wrapper(mobile, reference):
        """
        Static method to run the TM-Score command and return the TM-Score value.

        Parameters:
        ----------
        mobile : str
            The path to the mobile structure (the structure being compared).
        reference : str
            The path to the reference structure.

        Returns:
        -------
        tmscore : float
            The TM-Score value from the comparison.
        """
        command = f"TMscore {mobile} {reference} -outfmt 2"

        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        tmscore = None

        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                match = re.search(r'TM-score\s+=\s+([0-9.]+)', output)
                if match:
                    tmscore = float(match.group(1))
        process.poll()
        return tmscore

    def run_tmscore(self, folder_path, custom_ref):
        """
        Run TM-Score on all PDB files in a folder against a custom reference structure.

        Parameters:
        ----------
        folder_path : str
            The path to the folder containing PDB files.
        custom_ref : str
            The custom reference structure file path.

        Returns:
        -------
        tmscore_dict : dict
            A dictionary containing the TM-Score results for each frame.
        """
        pdb_files = sorted(glob(os.path.join(folder_path, '*.pdb')))
        if not custom_ref:
            custom_ref = pdb_files[0]
        tmscore_dict = {'frame': [], 'tmscore': []}
        for idx, pdb_file in enumerate(pdb_files):
            tmscore_dict['frame'].append(idx)
            tmscore_dict['tmscore'].append(self.tmscore_wrapper(pdb_file, custom_ref))
        return tmscore_dict

    def calculate_2d_tmscore(self, trial):
        """
        Calculate 2D TM-Score for a given trial.

        Parameters:
        ----------
        trial : str
            The identifier for the trial being analyzed.

        Returns:
        -------
        tmscore_2d_data : np.ndarray
            A 2D array of TM-Score values against two reference structures.
        """
        predictions_path_trial = f"{self.input_dict['predictions_path']}/{trial}"
        universe = self.prediction_dicts[trial]['mda_universe']

        ref_paths = (self.ref_gr, self.ref_alt)

        if self.slice_predictions:
            predictions_path_trial = f"{self.input_dict['predictions_path']}/tmp_pdb"
            create_directory(predictions_path_trial)
            self.slice_models(universe, self.slice_predictions, predictions_path_trial)

            ref_paths = []
            for idx, reference in enumerate((self.ref_gr, self.ref_alt)):
                universe = mda.Universe(reference)
                ag = universe.select_atoms(self.slice_predictions)
                ref_path = f"{predictions_path_trial}/ref_{idx}.pdb"
                ref_paths.append(ref_path)
                ag.write(ref_path)

        tmscore_gr = self.run_tmscore(predictions_path_trial, ref_paths[0])
        tmscore_alt = self.run_tmscore(predictions_path_trial, ref_paths[1])
        tmscore_2d_data = np.array([tmscore_gr['tmscore'], tmscore_alt['tmscore']])

        if self.slice_predictions:
            shutil.rmtree(predictions_path_trial)
        return tmscore_2d_data

    def fit_and_filter_data(self, tmscore_2d_data, n_stdevs):
        """
        Fit a parabola to the 2D TM-Score data and filter points based on the standard deviation threshold.

        Parameters:
        ----------
        tmscore_2d_data : np.ndarray
            A 2D array of TM-Score values.
        n_stdevs : int
            Number of standard deviations to use for filtering the data.

        Returns:
        -------
        None
        """
        popt, _ = curve_fit(parabola, tmscore_2d_data[0], tmscore_2d_data[1])
        fit_x = np.linspace(min(tmscore_2d_data[0]), max(tmscore_2d_data[0]), 100)
        fit_y = parabola(fit_x, *popt)

        fitted_curve_values = parabola(tmscore_2d_data[0], *popt)
        distances = np.abs(tmscore_2d_data[1] - fitted_curve_values)

        threshold = np.mean(distances) + n_stdevs * np.std(distances)
        close_points = distances < threshold
        x_close = tmscore_2d_data[0, close_points]
        y_close = tmscore_2d_data[1, close_points]

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

    def plot_filtering_data(self, tmscore_2d_data):
        """
        Generate and save a plot of the filtered 2D TM-Score data with the fitted curve.

        Parameters:
        ----------
        tmscore_2d_data : np.ndarray
            A 2D array of TM-Score values.

        Returns:
        -------
        None
        """
        plt.figure(figsize=(5, 4))
        plt.scatter(tmscore_2d_data[0], tmscore_2d_data[1], s=10)
        plt.plot(self.filtering_dict['fit_x'], self.filtering_dict['fit_y'], label='Fitted Curve', color='red')
        plt.scatter(self.filtering_dict['x_close'], self.filtering_dict['y_close'], label='Close Points', color='green', s=20)

        title = f"{self.input_dict['jobname']} {self.input_dict['max_seq']} {self.input_dict['extra_seq']}"
        plt.title(title, fontsize=15)
        plt.xlabel('TM-Score vs. Ref1', fontsize=14)
        plt.ylabel('TM-Score vs. Ref2', fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.legend(loc='best')
        plt.tight_layout()

        plot_path = (f"{self.input_dict['output_path']}/"
                     f"{self.input_dict['jobname']}/"
                     f"analysis/"
                     f"tmscore_2d/"
                     f"{self.input_dict['jobname']}_"
                     f"{self.input_dict['max_seq']}_"
                     f"{self.input_dict['extra_seq']}_"
                     f"tmscore_2d_fit.png")

        plt.savefig(plot_path, dpi=300)
        plt.close()

    def cluster_2d_data(self, tmscore_2d_data, n_clusters):
        """
        Perform clustering on the filtered 2D TM-Score data and store clustering results.

        Parameters:
        ----------
        tmscore_2d_data : np.ndarray
            A 2D array of TM-Score values.
        n_clusters : int
            Number of clusters to form.

        Returns:
        -------
        None
        """
        kmeans = KMeans(n_clusters)
        close_points_2d = np.array([self.filtering_dict['x_close'], self.filtering_dict['y_close']]).T

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
        outliers = total_samples - tmscore_2d_data.shape[0]
        outliers = 100 - ((outliers / total_samples) * 100)

        for k in unique_labels:
            class_member_mask = (correct_labels == k)
            xy = close_points_2d[class_member_mask]
            cluster_counts[k] = round(((len(xy)) / total_samples) * 100, 1)

        self.clustering_dict = {
            'labels': labels,
            'correct_labels': correct_labels,
            'centroids': centroids,
            'cluster_counts': cluster_counts,
            'outliers': outliers,
            'unique_labels': unique_labels,
            'close_points_2d': close_points_2d
        }

    def plot_and_save_2d_data(self, output_path):
        """
        Plot the clustered 2D TM-Score data, save the plot, and return a DataFrame with the clustering information.

        Parameters:
        ----------
        output_path : str
            The path where the plot will be saved.

        Returns:
        -------
        df : pd.DataFrame
            A DataFrame containing the clustering information.
        """
        colors = ['blue', 'green', 'magenta', 'orange', 'grey', 'brown', 'cyan', 'purple']
        unique_labels = self.clustering_dict['unique_labels']
        correct_labels = self.clustering_dict['correct_labels']
        cluster_counts = self.clustering_dict['cluster_counts']
        centroids = self.clustering_dict['centroids']
        outliers = self.clustering_dict['outliers']

        plotter = self.widget.add_plot(centroids[:, 0], centroids[:, 1], title=title, x_label='TM-Score vs. Ref1 (Å)', y_label='TM-Score vs. Ref2 (Å)', label='Centroids', scatter=True)
        self.widget.add_line(plotter, self.filtering_dict['fit_x'], self.filtering_dict['fit_y'], label='Fitted Curve', color='r')
    
        for i in unique_labels:
            cluster_points = self.clustering_dict['close_points_2d'][correct_labels == i]
            self.widget.add_scatter(plotter, cluster_points[:, 0], cluster_points[:, 1], color=colors[i], label=f'Cluster {i} pop: {cluster_counts[i]}')
  
        if output_path:
            plt.figure(figsize=(5, 4))
            for i in unique_labels:
                cluster_points = self.clustering_dict['close_points_2d'][correct_labels == i]
                plt.scatter(cluster_points[:, 0], cluster_points[:, 1], c=colors[i],
                            label=f'Cluster {i} pop: {cluster_counts[i]}', alpha=0.6)
            plt.scatter(centroids[:, 0], centroids[:, 1], s=100, c='black', marker='X', label='Centroids')
            plt.plot(self.filtering_dict['fit_x'], self.filtering_dict['fit_y'], label='Fitted Curve', color='red')

            title = (f"{self.input_dict['jobname']} {self.input_dict['max_seq']} {self.input_dict['extra_seq']} "
                    f"Score: {self.filtering_dict['ratio']:.2f}")

            plt.title(title, fontsize=16)
            plt.legend()
            plt.xlabel('TM-Score vs. Ref1', fontsize=14)
            plt.ylabel('TM-Score vs. Ref2', fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=12)
            plt.tight_layout()

            plot_path = (f"{self.input_dict['output_path']}/"
                        f"{self.input_dict['jobname']}/"
                        f"analysis/"
                        f"tmscore_2d/"
                        f"{self.input_dict['jobname']}_"
                        f"{self.input_dict['max_seq']}_"
                        f"{self.input_dict['extra_seq']}_"
                        f"tmscore_2d_clustered.png")

            plt.savefig(plot_path, dpi=300)
            plt.close()

        records = []
        for i in unique_labels:
            records.append({
                'trial': self.input_dict['trial'],
                'score': self.filtering_dict['ratio'],
                'cluster_label': i,
                'cluster_pop': cluster_counts[i],
                'centroid_values': centroids[i],
                '%_outliers': outliers
            })

        df = pd.DataFrame(records)
        return df

    def get_2d_tmscore(self, tmscore_mode_df_path, n_stdevs, n_clusters, output_path):
        """
        Execute the full 2D TM-Score analysis for all trials, including fitting, filtering, clustering, and saving results.

        Parameters:
        ----------
        tmscore_mode_df_path : str
            The path to the TM-Score mode data file.
        n_stdevs : int
            Number of standard deviations to use for filtering the data.
        n_clusters : int
            Number of clusters to form.
        output_path : str
            The path where the results will be saved.

        Returns:
        -------
        None
        """
        df = pd.read_csv(tmscore_mode_df_path)
        unique_trials = df['trial'].unique()

        df_all_trials = pd.DataFrame()
        with tqdm(total=len(self.prediction_dicts), bar_format=TQDM_BAR_FORMAT) as pbar:
            for trial in unique_trials:
                if not n_clusters:
                    unique_df = df[df['trial'] == trial]
                    n_clusters_trial = len(unique_df['mode_label']) + 1
                pbar.set_description(f'Running 2D TM-Score analysis for {trial}')
                self.input_dict['trial'] = trial
                self.input_dict['max_seq'] = self.prediction_dicts[trial]['max_seq']
                self.input_dict['extra_seq'] = self.prediction_dicts[trial]['extra_seq']

                tmscore_2d_data = self.calculate_2d_tmscore(trial)
                if len(tmscore_2d_data) > 0:
                    self.fit_and_filter_data(tmscore_2d_data, n_stdevs)
                    self.plot_filtering_data(tmscore_2d_data)
                    self.show_filt_data(self, tmscore_2d_data)
                    self.cluster_2d_data(tmscore_2d_data, n_clusters_trial)
                    df_to_save = self.plot_and_save_2d_data(output_path)
                    df_all_trials = pd.concat([df_all_trials, df_to_save], ignore_index=True)

                pbar.update(n=1)

            csv_path = (f"{self.input_dict['output_path']}/"
                        f"{self.input_dict['jobname']}/"
                        f"analysis/"
                        f"tmscore_2d/"
                        f"{self.input_dict['jobname']}_"
                        f"clustering_analysis.csv")


        print(f"\nSaving {self.input_dict['jobname']} 2D TM-Score analysis results to {csv_path}\n")

        df_all_trials.to_csv(csv_path, index=False)
