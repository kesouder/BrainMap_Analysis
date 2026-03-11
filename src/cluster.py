from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_hierarchical_clustering(corr_matrix, labels, method="average", figsize=(12, 8), n_cluster=2):
    """
    Performs hierarchical clustering on a correlation matrix and plots the dendrogram.

    Parameters:
        corr_matrix (DataFrame or 2D array-like): Correlation matrix for clustering.
        labels (list): Labels for the dendrogram leaves.
        method (str): Linkage method for clustering (e.g., "average", "complete", "ward").
        figsize (tuple): Size of the dendrogram figure.
    """
    # Distance matrix from correlation
    D = 1 - corr_matrix.values
    np.fill_diagonal(D, 0)

    # Condensed distance for linkage
    D_condensed = squareform(D, checks=False)

    # Perform hierarchical clustering
    Z = linkage(D_condensed, method=method)

    # Plot dendrogram
    plt.figure(figsize=figsize)
    plt.title("Dendrogram of Hierarchical Clustering")
    dendrogram(Z, labels=labels, leaf_rotation=45)
    plt.tight_layout()
    plt.savefig('images/Dendrogram.jpg')
    plt.show()
    n_clusters = n_cluster
    h_labels = fcluster(Z, t=n_clusters, criterion="maxclust")
    return Z, h_labels

def get_cluster_df(corr_df, h_labels):
    h_cluster_df = pd.DataFrame({"map": corr_df.index, "cluster": h_labels}).set_index("map")
    h_cluster_df.sort_values("cluster")
    return h_cluster_df