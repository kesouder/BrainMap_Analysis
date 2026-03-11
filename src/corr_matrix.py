import pandas as pd
import numpy as np
import seaborn as sns
from neuromaps import nulls, images, transforms, datasets, resampling
from scipy.stats import zscore, rankdata
import matplotlib.pyplot as plt
def is_volumetric(map_dict):
    """Check if a map is volumetric (MNI space)"""
    return map_dict.get('space') in ['MNI152', 'MNI305']

## Load and transform the dara
def load_and_prepare_map(map_dict, readable_map_names=None, target_space='fsLR', target_den='32k'):
    """
    Load a brain map and prepare it for analysis
    - Volumetric maps: transform to surface
    - Surface maps: resample if needed, use both hemispheres
    - return (transformed data array, space, density)
    """
    # Determine display name for logging
    desc = map_dict['desc']
    if readable_map_names and desc in readable_map_names:
        display_name = readable_map_names[desc]
    else:
        display_name = desc
    
    # Fetch the map
    brain_map = datasets.fetch_annotation(**map_dict)
    
    # Check if volumetric
    if is_volumetric(map_dict):
        print(f"Volumetric map detected, transforming to {target_space} surface...")
        
        # Transform volumetric to surface (both hemispheres)
        try:
            surface_map = transforms.mni152_to_fslr(brain_map, fslr_density=target_den)
            # surface_map will be a tuple of (left_hemi, right_hemi)
            lh_data = images.load_data(surface_map[0])
            rh_data = images.load_data(surface_map[1])
            full_data = np.hstack([lh_data, rh_data])
            return full_data, target_space, target_den
            
        except Exception as e:
            print(f"    Warning: Could not transform volumetric map: {e}")
 
    # Surface map processing
    else:
        src_space = map_dict.get('space')
        src_den = map_dict.get('den')
        template = {'source': 'hcps1200', 'desc': 'thickness', 'space': 'fsLR', 'den': '32k'}
        template_map = datasets.fetch_annotation(**template)
        # Load both hemispheres
        if len(brain_map) == 2:
            lh_map, rh_map = brain_map
        else:
            raise ValueError(f"Expected 2 hemispheres for surface map, got {len(brain_map)}")
        
        if src_den != target_den or src_space != target_space:
            print(f"Resampling from {src_space}-{src_den} to {target_space}-{target_den}...")            
            # Resample both hemispheres
            lh_resampled, _ = resampling.resample_images(
                lh_map, template_map[0],
                src_space=src_space, trg_space=target_space,
                hemi='L', resampling='transform_to_trg'
            )
            rh_resampled, _ = resampling.resample_images(
                rh_map, template_map[1],
                src_space=src_space, trg_space=target_space,
                hemi='R', resampling='transform_to_trg'
            )
            
            lh_data = images.load_data(lh_resampled)
            rh_data = images.load_data(rh_resampled)
        else:
            lh_data = images.load_data(lh_map)
            rh_data = images.load_data(rh_map)
        
        full_data = np.hstack([lh_data, rh_data])
        return full_data, src_space if src_den == target_den else target_space, \
               src_den if src_den == target_den else target_den
    
def compute_spin_spearman_analysis(
                                    brain_maps,
                                    n_perm=1000,
                                    target_space='fsLR',
                                    target_den='32k',
                                    seed=1234,
                                    readable_map_names=None):
    """
    Spearman-only brain map analysis with Alexander-Bloch spin correction.

    Returns:
        cov_empirical      (Spearman-based covariance matrix)
        cov_pvalues        (Spin-based p-values)
        corr_empirical     (Spearman correlation matrix)
        df_corr
        df_pvals
        df_significant
        map_names
        prepared_maps
        nulls_list
    """
    prepared_maps = []
    map_names = []

    for brain_map in brain_maps:
        prep_map, final_space, final_den = load_and_prepare_map(
            brain_map, readable_map_names
        )
        prepared_maps.append(prep_map)
        map_names.append(f"{brain_map['source']}_{brain_map['desc']}")

    n_maps = len(prepared_maps)

    nulls_list = []
    for prep_map in prepared_maps:
        rotated = nulls.alexander_bloch(
            prep_map,
            atlas=target_space,
            density=target_den,
            n_perm=n_perm,
            seed=seed
        )
        nulls_list.append(rotated)


    data_stack = np.column_stack(prepared_maps)
    global_mask = np.all(~np.isnan(data_stack), axis=1)

    clean_data = data_stack[global_mask]
    n_vertices_clean = clean_data.shape[0]


    ranked = np.apply_along_axis(rankdata, 0, clean_data)
    z_ranked = zscore(ranked, axis=0)

    corr_empirical = (z_ranked.T @ z_ranked) / n_vertices_clean

    p_counts = np.zeros((n_maps, n_maps))
    clean_nulls_list = [n[global_mask, :] for n in nulls_list]

    for k in range(n_perm):
        null_data_k = np.zeros_like(z_ranked)

        for m in range(n_maps):
            spin_k = clean_nulls_list[m][:, k]
            spin_ranked = rankdata(spin_k)
            null_data_k[:, m] = zscore(spin_ranked)

        null_corr_k = np.corrcoef(null_data_k.T)

        p_counts += (np.abs(null_corr_k) >= np.abs(corr_empirical))

    cov_pvalues = (p_counts + 1) / (n_perm + 1)

    stds = np.std(clean_data, axis=0)
    cov_empirical = corr_empirical * np.outer(stds, stds)


    display_names = map_names
    if readable_map_names is not None:
        display_names = [
            readable_map_names.get(name, name) for name in map_names
        ]

    df_corr = pd.DataFrame(corr_empirical, index=display_names, columns=display_names).round(3)
    df_pvals = pd.DataFrame(cov_pvalues, index=display_names, columns=display_names).round(3)

    df_significant = df_corr.copy()
    df_significant[df_pvals > 0.05] = 0

    return (
        cov_empirical,
        cov_pvalues,
        corr_empirical,
        df_corr,
        df_pvals,
        df_significant,
        map_names,
        prepared_maps,
        nulls_list
    )


## Function to plot heatmaps
def plot_correlation_heatmaps(df_corr, df_significant, method=None):
    """
    Plots two side-by-side heatmaps: 
    1. The Full Correlation Matrix
    2. The Significant-Only Matrix (Non-significant cells masked)
    """
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    
    sns.heatmap(
        df_corr, 
        annot=True, 
        fmt=".2f", 
        cmap='RdBu', 
        center=0, 
        vmin=-1, vmax=1,
        square=True,
        ax=axes[0],
        cbar_kws={'label': method}
    )
    axes[0].set_title("Full Correlation Matrix", fontsize=16)
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].tick_params(axis='y', rotation=0)
    
    # plot matrix masked
    mask = (df_significant == 0)
    
    sns.heatmap(
        df_significant, 
        annot=True, 
        fmt=".2f", 
        cmap='RdBu_r', 
        center=0, 
        vmin=-1, vmax=1,
        square=True,
        mask=mask, # hide non significant cells
        ax=axes[1],
        cbar_kws={'label': method}
    )
    axes[1].set_title("Significant Correlations Only (p < 0.05)", fontsize=16)
    axes[1].tick_params(axis='x', rotation=45)
    axes[1].tick_params(axis='y', rotation=0)
    
    plt.tight_layout()
    plt.show()

# Plot 1 Correlation Heatmap, no use of subplots
def plot_correlation_heatmap(corr_matrix, labels, title="Correlation Matrix Heatmap", figsize=(7, 6), cmap='RdBu'):
    """
    Plots a heatmap for a given correlation matrix.

    Parameters:
        corr_matrix (2D array-like): Correlation matrix to visualize.
        labels (list): Labels for the x and y axes.
        title (str): Title of the heatmap.
        figsize (tuple): Size of the figure.
        cmap (str): Colormap for the heatmap.
    """
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(corr_matrix, vmin=-1, vmax=1, cmap=cmap)  # correlation range

    # ticks + labels
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)

    # grid lines to separate cells
    ax.set_xticks(np.arange(-.5, len(labels), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
    ax.grid(which="minor", linestyle="-", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    # add colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Correlation")

    ax.set_title(title)
    plt.tight_layout()
    plt.show()
