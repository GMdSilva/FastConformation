import os
import pandas as pd
import numpy as np

from glob import glob

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.lines import Line2D
import seaborn as sns

from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from sklearn.cluster import KMeans

from MDAnalysis.analysis import align, pca

from analysis_utils import parabola

def pca_from_ensemble(jobname,
                      prediction_dicts,
                      output_path,
                      align_range,
                      analysis_ranges):

    pcas = {}
    n_clusters = 3

    for result in prediction_dicts:
        u = prediction_dicts[result]['mda_universe']
        max_seq = prediction_dicts[result]['max_seq']
        extra_seq = prediction_dicts[result]['extra_seq']

        for analysis_range_idx, analysis_range in enumerate(analysis_ranges):
            aligner = align.AlignTraj(u, u, select=align_range, in_memory=True).run()

            pc = pca.PCA(u, select=analysis_range,
                         align=False, mean=None,
                         n_components=None).run()

            backbone = u.select_atoms(analysis_range)

            # Transform and select only the first two components
            transformed = pc.transform(backbone, n_components=3)
            df = pd.DataFrame(transformed, columns=['PC1', 'PC2', 'PC3'])

            kmeans = KMeans(n_clusters=n_clusters)
            df['Cluster'] = kmeans.fit_predict(df[['PC1', 'PC2']])
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
            colors = ['blue', 'purple', 'green']

            sns.scatterplot(x='PC1', y='PC2', hue='Cluster', palette='viridis', data=df, legend=False)
            plt.plot(df['PC1'], df['PC2_fit'], color='red')  # Parabola fit line
            plt.title(f"{jobname} Parabola Fit R²: {r2:.2f}", fontsize=18)
            plt.tick_params(axis='both', which='major', labelsize=15)
            plt.xlabel('PC1', fontsize=17)
            plt.ylabel('PC2', fontsize=17)
            cluster_populations = df['Cluster'].value_counts().sort_index()
            legend_handles = [Line2D([0], [0], marker='o', color='w',
                                     label=f'Cluster {j} (Pop: {cluster_populations.get(j, 0)})',
                                     markersize=10, markerfacecolor=colors[j]) for j in range(n_clusters)]

            plt.legend(handles=legend_handles, fontsize=10, loc='upper right')

            plt.show()

    return pcas


def order_frames_by_pca(base_path, prefix, align_range, analysis_range, sort=True):
    pcas = {}
    n_clusters = 2

    for dir in os.listdir(base_path):
        if os.path.isdir(os.path.join(base_path, dir)) and dir.startswith(prefix):
            aligner = align.AlignTraj(u, u, select=align_range,
                                      in_memory=True).run()

            pc = pca.PCA(u, select=analysis_range,
                         align=False, mean=None,
                         n_components=None).run()

            backbone = u.select_atoms(analysis_range)

            # Transform and select only the first two components
            transformed = pc.transform(backbone, n_components=3)
            order = np.argsort(transformed[:, 0])
            #file_list = glob.glob('kit_results/predictions/kinase/5121024/test/kit_d816h_5121024/*.pdb')
            file_list = glob.glob('abl_results/predictions/kinase/5121024/abl_wt_5121024/*.pdb')
            print(file_list)
            ordered_file_list = [file_list[i] for i in order]
            topology_file = pdb_files[0]

            u_getter = loadPredictions(ordered_file_list, topology_file)
            u_getter.load()
            u = u_getter.get_universe()
            #ag = u.select_atoms('(backbone and resid 1-88) or (backbone and resid 178-346)')
            ag = u.select_atoms('backbone')
            aligner = align.AlignTraj(u, u, select=align_range,
                                      in_memory=True).run()
            ag.write("abl_c-alpha.pdb", frames='all')
            return order


def pca_from_ensemble_movie(jobname,
                      prediction_dicts,
                      output_path,
                      align_range,
                      analysis_ranges):

    pcas = {}
    n_clusters = 3

    for result in prediction_dicts:
        u = prediction_dicts[result]['mda_universe']
        max_seq = prediction_dicts[result]['max_seq']
        extra_seq = prediction_dicts[result]['extra_seq']

        for analysis_range_idx, analysis_range in enumerate(analysis_ranges):
            aligner = align.AlignTraj(u, u, select=align_range, in_memory=True).run()

            pc = pca.PCA(u, select=analysis_range,
                         align=False, mean=None,
                         n_components=None).run()

            backbone = u.select_atoms(analysis_range)

            # Transform and select only the first two components
            transformed = pc.transform(backbone, n_components=3)
            df = pd.DataFrame(transformed, columns=['PC1', 'PC2', 'PC3'])
            print(df)

            kmeans = KMeans(n_clusters=n_clusters)
            df['Cluster'] = kmeans.fit_predict(df[['PC1', 'PC2']])
            pcas[f'{result}_{analysis_range}'] = df
            order = np.argsort(transformed[:, 0])
            df = df.iloc[order]

            # Set up the plot
            fig, ax = plt.subplots()
            colors = sns.color_palette("viridis", n_clusters)

            # Determine fixed axis ranges
            x_range = [df['PC1'].min(), df['PC1'].max()]
            y_range = [df['PC2'].min(), df['PC2'].max()]

            def animate(i):
                if i == 0:
                    ax.clear()
                    ax.set_xlim(x_range)
                    ax.set_ylim(y_range)
                    ax.tick_params(axis='both', which='major', labelsize=13)
                    ax.set_xlabel('PC1', fontsize=15)
                    ax.set_ylabel('PC2', fontsize=15)

                ax.set_title(f"{jobname} Prediction - Frame {i}", fontsize=17)

                if i >= 1:
                    point = df.iloc[i - 1]
                    cluster = int(point['Cluster'])  # Cast cluster to integer
                    ax.scatter(point['PC1'], point['PC2'], c=[colors[cluster]])

                # Calculate and update cluster populations up to the current frame
                cluster_populations = df.iloc[:i]['Cluster'].value_counts().sort_index()

                # Create custom legend handles
                legend_handles = [Line2D([0], [0], marker='o', color='w',
                                         label=f'Cluster {j} (Pop: {cluster_populations.get(j, 0)})',
                                         markersize=10, markerfacecolor=colors[j]) for j in range(n_clusters)]

                ax.legend(handles=legend_handles, fontsize=10, loc='upper right')

            anim = FuncAnimation(fig, animate, frames=len(df) + 1, interval=10, repeat=False)

            # Save the animation
            anim.save(os.path.join(output_path,
                                   f'{jobname}_{max_seq}_{extra_seq}_{analysis_range_idx}_pca_animation.gif'),
                      writer='Pillow')

            return pca