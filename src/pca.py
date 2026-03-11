import numpy as np
import matplotlib.pyplot as plt

def plot_pca_and_scree(corr_matrix, labels, figsize_scatter=(6, 6), figsize_scree=(6, 4)):
    """
    Performs PCA on a correlation matrix and plots the loadings scatter plot and scree plot.

    Parameters:
        corr_matrix (DataFrame or 2D array-like): Correlation matrix for PCA.
        labels (list): Labels for the scatter plot points.
        figsize_scatter (tuple): Size of the scatter plot figure.
        figsize_scree (tuple): Size of the scree plot figure.
    """
    # Perform PCA
    eigvals, eigvecs = np.linalg.eigh(corr_matrix.values)  # ascending
    idx = np.argsort(eigvals)[::-1]  # sort descending
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    
    # Variance explained
    var_exp = eigvals / eigvals.sum()
    
    # Calculate Loadings
    loadings = eigvecs * np.sqrt(eigvals)
    
    # Coordinates for PC1/PC2 plot
    x = loadings[:, 0]
    y = loadings[:, 1]
    
    # Plot: loadings scatter with labels
    plt.figure(figsize=figsize_scatter)
    plt.axhline(0, linewidth=1)
    plt.axvline(0, linewidth=1)
    plt.scatter(x, y)
    for i, name in enumerate(labels):
        plt.text(x[i], y[i], "  " + name, va="center")
    plt.xlabel(f"PC1 loading (var={var_exp[0]:.2%})")
    plt.ylabel(f"PC2 loading (var={var_exp[1]:.2%})")
    plt.title("PCA of Correlation Matrix (Loadings)")
    plt.tight_layout()
    plt.show()
    
    # Plot: Scree plot
    plt.figure(figsize=figsize_scree)
    plt.plot(np.arange(1, len(var_exp) + 1), var_exp, marker="o")
    plt.xticks(np.arange(1, len(var_exp) + 1))
    plt.xlabel("Principal component")
    plt.ylabel("Proportion variance explained")
    plt.title("Scree Plot (from Correlation Matrix)")
    plt.tight_layout()
    plt.show()
    return loadings, var_exp

def plot_pca_loadings_heatmap(loadings, labels, var_exp, figsize=(10, 6), cmap='RdBu'):
    """
    Plots a heatmap of PCA loadings.

    Parameters:
        loadings (2D array-like): PCA loadings matrix (PCs × features).
        labels (list): Labels for the columns (features).
        var_exp (array-like): Variance explained by each principal component.
        figsize (tuple): Size of the heatmap figure.
        cmap (str): Colormap for the heatmap.
    """
    L = loadings.T  # Transpose to (PCs × features)
    plt.figure(figsize=figsize)
    im = plt.imshow(L, aspect="auto", vmin=-1, vmax=1, cmap=cmap, alpha=0.8)
    
    # Axis labels
    plt.xticks(np.arange(len(labels)), labels, rotation=45, ha="right")
    plt.yticks(np.arange(len(L)), [f"PC{i+1} ({var_exp[i]:.1%})" for i in range(len(L))])
    plt.title("PCA Loadings Heatmap (PCs × Features)")
    
    # Annotate each cell value
    for i in range(len(L)):
        for j in range(len(labels)):
            plt.text(j, i, f"{L[i, j]:.2f}", ha="center", va="center", fontsize=8)
    
    plt.colorbar(im, label="Loading")
    plt.tight_layout()
    plt.show()
# L = loadings.T  # (8 PCs, 8 maps)
# plt.figure(figsize=(10, 6))
# im = plt.imshow(L, aspect="auto", vmin=-1, vmax=1, cmap='RdBu', alpha=0.8)
# # axis labels
# plt.xticks(np.arange(len(map_names)), map_names, rotation=45, ha="right")
# plt.yticks(np.arange(len(L)), [f"PC{i+1} ({var_exp[i]:.1%})" for i in range(len(L))])
# plt.title("PCA loadings heatmap (PCs × maps)")
# # annotate each cell value
# for i in range(len(L)):
#     for j in range(len(map_names)):
#         plt.text(j, i, f"{L[i, j]:.2f}", ha="center", va="center", fontsize=8)
# plt.colorbar(im, label="Loading")
# plt.tight_layout()
# plt.show()