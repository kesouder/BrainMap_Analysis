import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, rankdata, spearmanr
import math
## Funciton to plot EDA pairwise scatter plots 

def plot_pairwise_scatterplots(prepared_maps, map_names, readable_map_names=None):
    # data prep
    # 32k rows 8 columns
    data_stack = np.column_stack(prepared_maps)
    global_mask = np.all((~np.isnan(data_stack)) & (data_stack != 0), axis=1) # get rid of NaNs
    clean_data = data_stack[global_mask]
    
    n_maps = len(prepared_maps)
    
    # making sure brain maps names exist and are readable
    # uses readable names instead of desc name
    display_names = map_names
    if readable_map_names:
        display_names = []
        for name in map_names:
            if name in readable_map_names:
                display_names.append(readable_map_names[name])
                continue
            found_match = False
            for key, val in readable_map_names.items():
                if name.endswith(f"_{key}"):
                    display_names.append(val)
                    found_match = True
                    break
            if not found_match:
                display_names.append(name)

    # 28 unique pairs
    n_pairs = (n_maps * (n_maps - 1)) // 2
    n_cols = 4
    n_rows = math.ceil(n_pairs / n_cols)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 4 * n_rows))
    axes = axes.flatten()
    
    plot_idx = 0
    
    for i in range(n_maps):
        for j in range(i + 1, n_maps):
            ax = axes[plot_idx]
            
            x = clean_data[:, i] # column i for all j rows
            y = clean_data[:, j]
            r, _ = pearsonr(x, y)
            rho, _ = spearmanr(x, y)
            # use hexbin for density instead of scatter plot because of ~32k points
            hb = ax.hexbin(x, y, gridsize=50, cmap='Blues', mincnt=1, bins='log')
            # add trend line
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            x_range = np.linspace(x.min(), x.max(), 100)
            ax.plot(x_range, p(x_range), "r--", alpha=0.8, linewidth=2)
            
            ax.set_xlabel(display_names[i], fontsize=8)
            ax.set_ylabel(display_names[j], fontsize=8)
            
            # Green if metrics agree (diff < 0.1), Red if they diverge
            diff = abs(r - rho)
            title_color = 'green' if diff < 0.1 else 'red'
            
            ax.set_title(f"r={r:.2f} | ρ={rho:.2f}", fontsize=10, color=title_color, fontweight='bold')
            ax.set_xticks([])
            ax.set_yticks([])
            
            plot_idx += 1
        
    for k in range(plot_idx, len(axes)):
        axes[k].axis('off')
        
    plt.tight_layout()
    plt.savefig("images/pairwise_scatter_plot_unranked.jpg")
    plt.show()

def pairwise_rank_scatterplot(prepared_maps, map_names, readable_map_names=None):
    data_stack = np.column_stack(prepared_maps)
    global_mask = np.all((~np.isnan(data_stack)) & (data_stack != 0), axis=1)
    clean_data = data_stack[global_mask]
    
    # Rank Transform - converts raw values to ranks (0, 1, 2... N)
    ranked_data = np.apply_along_axis(rankdata, 0, clean_data)
    
    n_maps = len(prepared_maps)
    display_names = map_names
    if readable_map_names:
        display_names = []
        for name in map_names:
            if name in readable_map_names:
                display_names.append(readable_map_names[name])
                continue
            found_match = False
            for key, val in readable_map_names.items():
                if name.endswith(f"_{key}"):
                    display_names.append(val)
                    found_match = True
                    break
            if not found_match: display_names.append(name)

    n_pairs = (n_maps * (n_maps - 1)) // 2
    n_cols = 4
    n_rows = math.ceil(n_pairs / n_cols)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 4 * n_rows))
    axes = axes.flatten()
    
    plot_idx = 0
    
    for i in range(n_maps):
        for j in range(i + 1, n_maps):
            ax = axes[plot_idx]
            x = ranked_data[:, i]
            y = ranked_data[:, j]
            
            # Compute stats (On ranks, Pearson == Spearman)
            rho, _ = pearsonr(x, y) 
            
            # Hexbin Plot
            hb = ax.hexbin(x, y, gridsize=50, cmap='viridis', mincnt=1, bins='log')
            
            # Trendline
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            x_range = np.linspace(x.min(), x.max(), 100)
            ax.plot(x_range, p(x_range), "r--", alpha=0.8, linewidth=2)
            
            ax.set_xlabel(f"{display_names[i]} (Rank)", fontsize=8)
            ax.set_ylabel(f"{display_names[j]} (Rank)", fontsize=8)
            
            ax.set_title(f"Spearman ρ={rho:.2f}", fontsize=10, fontweight='bold')
            ax.set_xticks([])
            ax.set_yticks([])
            
            plot_idx += 1
        
    for k in range(plot_idx, len(axes)):
        axes[k].axis('off')
        
    plt.tight_layout()
    plt.savefig("images/pairwise_scatterplot_ranked.jpg")
    plt.show()
