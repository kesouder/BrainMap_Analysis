import numpy as np
from statsmodels.stats.multitest import multipletests
import pandas as pd
from scipy.stats import zscore, rankdata
import matplotlib.pyplot as plt
import seaborn as sns
# Mutltiple Tests correction table
def generate_standard_comparison_table(df_corr, df_pvals, alpha=0.05):
    """
    Computes Bonferroni and FDR BH and FDR BY corrections
    """
    
    mask_tri = np.tril(np.ones(df_pvals.shape), k=-1).astype(bool)
    rows, cols = np.where(mask_tri)
    
    map1_names = df_corr.index[rows]
    map2_names = df_corr.columns[cols]
    
    r_values = df_corr.values[mask_tri]
    p_values_uncorrected = df_pvals.values[mask_tri]
    
    # can use p-val for appendix additional tables
    reject_bonf, pvals_bonf, _, _ = multipletests(p_values_uncorrected, alpha=alpha, method='bonferroni')
    reject_fdr_bh, pvals_fdr_bh, _, _ = multipletests(p_values_uncorrected, alpha=alpha, method='fdr_bh')
    reject_fdr_by, pvals_fdr_by, _, _ = multipletests(p_values_uncorrected, alpha=alpha, method='fdr_by')
    
    results = pd.DataFrame({
        'Map 1': map1_names,
        'Map 2': map2_names,
        'Spearman': r_values,
        'P-val (Uncorrected)': p_values_uncorrected,
        'Significant (Uncorrected)': p_values_uncorrected < alpha,
        'Significant (Bonf)': reject_bonf,
        'Significant (FDR-BH)': reject_fdr_bh,
        'Significant (FDR-BY)': reject_fdr_by
    })
    
    float_cols = ['Spearman', 'P-val (Uncorrected)']
    results[float_cols] = results[float_cols].round(4)
    
    return results

# Helper function to compute Max T threshold
def compute_maxt_threshold(prepared_maps, nulls_list, n_perm=1000, method='pearson'):
    data_stack = np.column_stack(prepared_maps)
    global_mask = np.all((~np.isnan(data_stack)) & (data_stack != 0), axis=1)
    
    clean_data = data_stack[global_mask]
    n_vertices_clean = clean_data.shape[0]
    clean_nulls_list = [n[global_mask, :] for n in nulls_list]
    
    max_null_stats = []

    for k in range(n_perm):        
        null_data_k = np.zeros((n_vertices_clean, len(prepared_maps)))
        for m in range(len(prepared_maps)):
            spin_k = clean_nulls_list[m][:, k]
            if method == 'spearman':
                spin_k = rankdata(spin_k, method='average')
            
            with np.errstate(divide='ignore', invalid='ignore'):
                z_spin = zscore(spin_k, nan_policy='omit')
            null_data_k[:, m] = z_spin

        null_data_k = np.nan_to_num(null_data_k)
        null_corr_matrix = (null_data_k.T @ null_data_k) / n_vertices_clean 
        np.fill_diagonal(null_corr_matrix, 0)
        
        max_r = np.nanmax(np.abs(null_corr_matrix))
        max_null_stats.append(max_r)
        
    threshold = np.nanpercentile(max_null_stats, 95)
    print(f"Max-T Critical Threshold: r > {threshold:.3f}")
    return threshold


# --- Main function for Max T table ---
def compute_and_display_maxt(df_corr, prepared_maps, nulls_list, method='pearson'):
    
    n_perm_actual = nulls_list[0].shape[1]
    
    # compute Threshold
    threshold = compute_maxt_threshold(prepared_maps, nulls_list, n_perm=n_perm_actual, method=method)
    
    if np.isnan(threshold):
        print("Error: Threshold is NaN.")
        return None, None

    mask_tri = np.tril(np.ones(df_corr.shape), k=-1).astype(bool)
    rows, cols = np.where(mask_tri)
    
    map1_names = df_corr.index[rows]
    map2_names = df_corr.columns[cols]
    r_values = df_corr.values[rows, cols]
    
    maxt_table = pd.DataFrame({
        'Map 1': map1_names,
        'Map 2': map2_names,
        'Correlation': r_values,
        'Significant': np.abs(r_values) > threshold  # Logic check applied to all
    }).round(3)
    
    maxt_table = maxt_table.sort_values(by='Correlation', key=abs, ascending=False).reset_index(drop=True)
        
    plt.figure(figsize=(10, 8))
    df_plot = df_corr.copy()
    visual_mask = np.abs(df_plot) <= threshold
    np.fill_diagonal(visual_mask.values, False)
    
    sns.heatmap(
        df_plot, 
        annot=True, 
        fmt=".2f", 
        cmap='RdBu_r', 
        center=0, 
        vmin=-1, vmax=1,
        square=True,
        mask=visual_mask,
        cbar_kws={'label': f'{method.capitalize()} Correlation'}
    )
    plt.title(f"Significant Couplings (Max-T Corrected)\nThreshold: |rho| > {threshold:.3f}", fontsize=14)
    plt.tight_layout()
    plt.savefig("images/max_t_corrected_table.jpg")
    plt.show()
    
    return maxt_table, threshold