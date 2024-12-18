import os

from sklearn.covariance import graphical_lasso,shrunk_covariance
from sklearn.preprocessing import StandardScaler
from sklearn.covariance import GraphicalLasso
from sklearn.covariance import EmpiricalCovariance
from sklearn.cluster import SpectralClustering, AgglomerativeClustering, DBSCAN
import numpy as np
import pandas as pd
import seaborn as sns
import umap
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import warnings
from sklearn.covariance import LedoitWolf, OAS
import sys

# Get the root directory
root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the Python directory to sys.path
python_dir = os.path.join(root_dir, 'Python')
sys.path.append(python_dir)

import gba_analysis

cluster_palette = ["#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", "#aa40fc", \
                  "#e377c2", "#b5bd61", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", \
                  "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31", \
                  "#b4d2b1", "#568f8b", "#1d4a60", "#cd7e59", "#ddb247", "#d15252", \
                  "#264653", "#2a9d8f", "#e9c46a", "#e76f51", "#f4a261", "#ef476f", \
                  "#ffd166","#06d6a0","#118ab2","#073b4c", "#fbf8cc","#fde4cf", \
                  "#ffcfd2","#f1c0e8","#cfbaf0","#a3c4f3","#90dbf4","#8eecf5", \
                  '#8359A3', '#5e503f', '#33CC99', '#F2C649', '#B94E48', '#0095B7', \
                  '#FF681F', '#e0aaff', '#FED85D', '#0a0908', '#C32148', '#98f5e1', \
                  "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", \
                  "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", \
                  "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", \
                  "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", \
                  "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F", \
                  "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", \
                  "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", \
                  "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C", \
                  "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81", \
                  "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", \
                  "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", \
                  "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", \
                  "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C", \
                  "#00B7FF", "#004DFF", "#00FFFF", "#826400", "#580041", "#FF00FF", "#00FF00", "#C500FF", \
                  "#B4FFD7", "#FFCA00", "#969600", "#B4A2FF", "#C20078", "#0000C1", "#FF8B00", "#FFC8FF", \
                  "#666666", "#FF0000", "#CCCCCC", "#009E8F", "#D7A870", "#8200FF", "#960000", "#BBFF00", \
                  "#FFFF00", "#006F00"]


path = './data/GBA2_v5/dsa_exp_download/'
diseases = {'DSA05072': 'Ulcerative colitis', 
    'DSA08946': 'Invasive Breast Carcinoma',
    'DSA01344': 'Asthma',
    'DSA08970': 'Lung Adenocarcinoma',
    'DSA08751': 'Dilated Cardiomyopathy',
    'DSA07529': 'Arthritis',
    'DSA07139': 'Psoriasis',
    'DSA09769': 'Chronic Obstructive Pulmonary Disease',
    'DSA00763': 'Type 2 diabetes'
}

# Specify the number of clusters
n_clusters = 40

for disease_id, disease in diseases.items():
    disease_name = disease.replace(' ', '_')
    profile_df = pd.read_csv(f'{path}{disease_id}_profile.csv', skiprows = 1)
    profile_df.set_index('GeneID', inplace=True)
    profile_df = profile_df.T

    # Scaling
    profile_df = pd.DataFrame(StandardScaler().fit_transform(profile_df), index=profile_df.index, columns=profile_df.columns)

    # Read metadata dataframe
    meta_df = pd.read_csv(f'{path}{disease_id}_meta.csv')
    meta_df.set_index('sample', inplace=True)
    meta_df.loc[profile_df.index, :]

    # Select the top highly variable genes
    gene_variances = profile_df.var()
    highly_variable = gene_variances.nlargest(2000).index
    highly_variable_df = profile_df[highly_variable]

    case_df = highly_variable_df.loc[meta_df['status'] == 'case'].T
    control_df = highly_variable_df.loc[meta_df['status'] == 'control'].T
    highly_variable_df = highly_variable_df.T
    print(case_df.shape, control_df.shape)

    # # Fit Graphical Lasso to estimate the precision matrix of case + control
    # glasso = GraphicalLasso(alpha=0.1, max_iter=200)  # Regularization parameter (tune as needed)
    # glasso.fit(highly_variable_df.T)
    # inv_cov_case_ctrl = glasso.precision_

    # # Fit Graphical Lasso to estimate the precision matrix of case + control
    # glasso = GraphicalLasso(alpha=0.1, max_iter=200)  # Regularization parameter (tune as needed)
    # glasso.fit(control_df.T)
    # inv_cov_ctrl = glasso.precision_

    # Step 3: Threshold the precision matrix to create an adjacency matrix
    #threshold = 0.1  # Set a threshold for sparsity
    #adjacency_matrix = (np.abs(precision_matrix) > threshold).astype(int)
    #np.fill_diagonal(adjacency_matrix, 0)  # Remove self-loops

    # Create graph
    # Disease + control
    # Fit Ledoit-Wolf covariance estimator
    lw = LedoitWolf()
    cov_case_ctrl = lw.fit(highly_variable_df.T).covariance_
    #cov_case_ctrl = np.cov(highly_variable_df, bias=False)
    try:
        inv_cov_case_ctrl = np.linalg.inv(cov_case_ctrl)
    except np.linalg.LinAlgError:
        print("Covariance matrix is singular, attempting regularized inversion")
        # Regularize with scikit-learn's EmpiricalCovariance if matrix is singular
        cov_estimator = EmpiricalCovariance().fit(cov_case_ctrl)
        inv_cov_case_ctrl = cov_estimator.precision_

    #  control
    # Fit Ledoit-Wolf covariance estimator
    lw = LedoitWolf()
    cov_ctrl = lw.fit(control_df.T).covariance_
    #cov_ctrl = np.cov(control_df, bias=False)
    try:
        inv_cov_ctrl = np.linalg.inv(cov_ctrl)
    except np.linalg.LinAlgError:
        print("Covariance matrix is singular, attempting regularized inversion")
        # Regularize with scikit-learn's EmpiricalCovariance if matrix is singular
        cov_estimator = EmpiricalCovariance().fit(cov_ctrl)
        inv_cov_ctrl = cov_estimator.precision_
        
    diff_df = inv_cov_case_ctrl - inv_cov_ctrl
    diff_df = diff_df.astype('float64')
    diff_df = np.abs(diff_df)
    diff_df = pd.DataFrame(np.maximum(diff_df, diff_df.T), index=highly_variable_df.index, columns=highly_variable_df.index ) # keep it symetric

    print(diff_df.shape)
    
    # Perform spectral clustering
    clustering = SpectralClustering(
        n_clusters=n_clusters,
        affinity='precomputed'  # Use 'precomputed' since we are providing an adjacency matrix
    #     assign_labels='kmeans'
    )

    # Fit the model and get cluster labels
    labels = clustering.fit_predict(diff_df)

    # UMAP can work with pre-computed distance matrices. Therefore, we calculate information hop distances from diff_df and paas it to umap.
    ihm_diff_df = gba_analysis.construct_information_hop_matrix(diff_df, diff_df.columns, verbose=True, use_max_distance=True)

    # The umap projection can not work if there are NaNs, +inf or -inf values or values that don't fit in float32,
    # so we need to perform the appropriate checks before applying umap and, if needed, perform re-scaling!    
    nans = ihm_diff_df.isna().sum()  # Check for NaNs
    plus_infs = (ihm_diff_df == float('inf')).sum()  # Check for positive infinity
    minus_infs = (ihm_diff_df == float('-inf')).sum()  # Check for negative infinity

    print(f"There are {nans} NaNs")
    print(f"There are {plus_infs} +inf values")
    print(f"There are {minus_infs} -inf values")

    # Check for values that are too large or too small to fit into float32
    float32_min = np.finfo(np.float32).min
    float32_max = np.finfo(np.float32).max

    # Check which values exceed the float32 range
    too_large = (ihm_diff_df > float32_max) | (ihm_diff_df < float32_min)
    problematic_value_count = too_large.sum().sum() # TODO: Fix this!

    # Print the columns and rows where values are too large
    print(f"There are {problematic_value_count} values in the matrix that don't fit into float32!")  # Shows how many problematic values in each column

    if (problematic_value_count > 0):
        print("Applying scaling for umap inputs.")
        # Initialize MinMaxScaler to scale data to a [0, 1] range
        scaler = MinMaxScaler()

        # Scale the matrix
        ihm_diff_df = scaler.fit_transform(ihm_diff_df)


    # Apply UMAP
    umap_model = umap.UMAP(n_components=2, random_state=42, metric="precomputed")  # For 2D projection
    umap_embedding = umap_model.fit_transform(ihm_diff_df)

    # Step 2: Create a DataFrame with UMAP coordinates and labels
    df_umap = pd.DataFrame(umap_embedding, columns=['UMAP1', 'UMAP2'])
    df_umap['label'] = labels

    # Step 3: Plot with Seaborn
    plt.figure(figsize=(10, 8))
    scatter_plot = sns.scatterplot(data=df_umap, x='UMAP1', y='UMAP2', hue='label', palette=cluster_palette[:n_clusters], s=50)
    plt.title(f'{disease}: UMAP Projection with cluster labels', size=30)
    plt.box(None)
    scatter_plot.legend(title='Cluster', bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)
    plt.savefig(f'{path}{disease_name}_exp_umap.png', dpi=200, bbox_inches='tight')
    plt.close()
    
    # Step 4: Plot UMAP per cluster
    nrows = 5
    fig, axs = plt.subplots(nrows, n_clusters//nrows, figsize=(30, 15))
    
    fig.suptitle(f'{disease}: Umap projection per cluster', size=32)
    xmax = df_umap['UMAP1'].max()
    xmin = df_umap['UMAP1'].min()
    ymax = df_umap['UMAP2'].max()
    ymin = df_umap['UMAP2'].min()
    label = 0
    clusters = {}
    conductances = []
    min_cond = 1
    min_cond_node = -1
    highlight_ax = None
    for row in range(nrows):
        for col in range(n_clusters//nrows):
            df_part = df_umap[labels==label]
            sns.scatterplot(data=df_part, x='UMAP1', y='UMAP2', color=cluster_palette[label], s=50, ax=axs[row, col])

            axs[row, col].set(yticklabels=[], xticklabels=[])
            axs[row, col].tick_params(bottom=False, left=False)
            axs[row, col].set_xlabel(None)
            axs[row, col].set_ylabel(None)
            axs[row, col].set_xlim(xmin,xmax)
            axs[row, col].set_ylim(ymin,ymax)

            clusters[label] = [highly_variable_df.index[i] for i, value in enumerate(labels) if value == label]
            cond = np.round(gba_analysis.calculate_conductance(clusters[label], (diff_df)), 3)
            conductances.append(cond)
            axs[row, col].set_title(f'{label} ({len(df_part)}): {cond}', size=20)
            if (cond < min_cond) and (len(df_part) <= 100) and (len(df_part) > 20):
                min_cond = cond
                min_cond_node = label
                highlight_ax = axs[row, col]  # Store the axis of the subplot to highlight

            label += 1
            #         sns.despine(bottom=True, left=True, ax=axs[row, col])
    
    # Highlight the subplot with the minimum conductance
    if highlight_ax:
        highlight_ax.patch.set_edgecolor('red')  # Set red border color
        highlight_ax.patch.set_linewidth(4)      # Set border thickness

    plt.savefig(f'{path}{disease_name}_exp_umap_per_cluster.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Smallest conductance in cluster: {min_cond_node} ({min_cond})')

    # Save graph and S
    diff_df = pd.DataFrame(diff_df, columns=case_df.index)
    g_path = f'{path}{disease_name}_exp_adj.csv'
    diff_df.to_csv(g_path, index=False)
    s_df = pd.DataFrame(clusters[min_cond_node], columns=['x'])
    s_path = f'{path}{disease_name}_exp_s.csv'
    s_df.to_csv(s_path, index=False)
    
    command = f"python ./tests/test_gba.py {disease_name} {g_path} {s_path} >> {path}gba2_{disease_name}.log"
    os.system(command)
    print("Done")
  