import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from sklearn.cluster import KMeans

from MDAnalysis.analysis import align, pca

from ensemble_analysis.analysis_utils import parabola

from tqdm import tqdm

TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'


def pca_from_ensemble(jobname,
                      prediction_dicts,
                      output_path,
                      align_range,
                      analysis_range,
                      n_clusters):

    pcas = {}
    print('')
    with tqdm(total=len(prediction_dicts), bar_format=TQDM_BAR_FORMAT) as pbar:
        for result in prediction_dicts:
            pbar.set_description(f'Running PCA for {result}')
            u = prediction_dicts[result]['mda_universe']
            max_seq = prediction_dicts[result]['max_seq']
            extra_seq = prediction_dicts[result]['extra_seq']

            align.AlignTraj(u, u, select=align_range, in_memory=True).run()

            pc = pca.PCA(u, select=analysis_range,
                         align=False, mean=None,
                         n_components=None).run()

            backbone = u.select_atoms(analysis_range)

            # Transform and select only the first two components
            transformed = pc.transform(backbone, n_components=3)
            df = pd.DataFrame(transformed, columns=['PC1', 'PC2', 'PC3'])

            kmeans = KMeans(n_clusters=n_clusters)
            df['Cluster'] = kmeans.fit_predict(df[['PC1', 'PC2']])

            labels = kmeans.labels_
            centroids = kmeans.cluster_centers_
            unique_labels = set(labels)
            cluster_counts = {i: 0 for i in unique_labels}
            for k in unique_labels:
                xy = df[df['Cluster'] == k]
                cluster_counts[k] = round(((len(xy)) / len(df['Cluster'])) * 100, 1)

            pcas[f'{result}_{analysis_range}'] = df
            order = np.argsort(transformed[:, 0])
            df = df.iloc[order]
            pcas[f'{result}_{analysis_range}'] = df

            # Curve fitting for PC1 and PC2 to a parabola
            popt, pcov = curve_fit(parabola, df['PC1'], df['PC2'])
            df['PC2_fit'] = parabola(df['PC1'], *popt)
            r2 = r2_score(df['PC2'], df['PC2_fit'])
            df.sort_values('PC1', inplace=True)

            # Create a scatter plot for PC1 and PC2 with clusters
            colors = ['blue', 'green', 'magenta', 'orange', 'grey', 'brown', 'cyan', 'purple']

            plt.figure(figsize=(5, 4))
            for i in unique_labels:
                cluster_points = df[df['Cluster'] == i]
                plt.scatter(cluster_points['PC1'], cluster_points['PC2'], c=colors[i],
                            label=f'Cluster {i} pop: {cluster_counts[i]}', alpha=0.6)
            plt.scatter(centroids[:, 0], centroids[:, 1], s=100, c='black', marker='X', label='Centroids')
            plt.plot(df['PC1'], df['PC2_fit'], label='Fitted Curve', color='red')

            plt.title(f"{jobname} {max_seq} {extra_seq} Fit R²: {r2:.2f}", fontsize=18)
            plt.tick_params(axis='both', which='major', labelsize=15)
            plt.xlabel('PC1', fontsize=17)
            plt.ylabel('PC2', fontsize=17)
            plt.legend(loc='best')

            full_output_path = f"{output_path}/{jobname}/analysis/pca/"
            figure_output_path = f"{full_output_path}/{jobname}_{max_seq}_{extra_seq}_pca.png"
            csv_output_path = f"{full_output_path}/{jobname}_{max_seq}_{extra_seq}_pca.csv"
            plt.tight_layout()
            plt.savefig(figure_output_path, dpi=300)
            df.to_csv(csv_output_path, index=False)
            pbar.update(n=1)

            plt.close()

    print(f"\nSaving {jobname} PCA analysis results to {full_output_path}\n")