import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr
from neuromaps.datasets import fetch_atlas
from neuromaps import plotting
from itertools import combinations
# Helper function
def get_readable_title(raw_name, name_dict):
    """Converts technical map name to readable title."""
    for key, readable in name_dict.items():
        if key in raw_name:
            return readable
    return raw_name

# Helper function
def generate_coupling_map(map1_data, map2_data, positive=True):
    """
    Generates a "Local Coupling Map"
    """
    mask = (~np.isnan(map1_data)) & (~np.isnan(map2_data)) & (map1_data != 0) & (map2_data != 0)
    coupling_map = np.zeros_like(map1_data, dtype=float)
    
    d1 = map1_data[mask]
    d2 = map2_data[mask]
    # Rank 0 to 1
    r1 = rankdata(d1, method='average') / len(d1)
    # check for negative relationship
    if positive:
        r2 = rankdata(d2, method='average') / len(d2)
    else:
        r2 = 1.0 - (rankdata(d2, method='average') / len(d2))

    diff = np.abs(r1 - r2)
    coupling_score = 1.0 - diff
    # cubed to account for noise
    coupling_score = coupling_score ** 3
    coupling_map[mask] = coupling_score
    return coupling_map

def plot_coupling_map(coupling_map_data, space='fsLR', den='32k', title='Local Coupling Strength', cmap='magma'):
    """
    Plots the coupling map using neuromaps plotting utilities.
    """
    atlas = fetch_atlas(space, den)
    surf_l = atlas['inflated'].L
    surf_r = atlas['inflated'].R
    n_vertices = len(coupling_map_data)
    mid_point = n_vertices // 2
    data_l = coupling_map_data[:mid_point]
    data_r = coupling_map_data[mid_point:]

    fig = plt.figure(figsize=(10, 4))
    vmin, vmax = 0, 1
    
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    plotting.plot_surf(
        surf_mesh=surf_l,
        surf_map=data_l,
        hemi='left',
        view='lateral',
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        axes=ax1,
        title='Left Hemisphere'
    )
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    plotting.plot_surf(
        surf_mesh=surf_r,
        surf_map=data_r,
        hemi='right',
        view='lateral',
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        axes=ax2,
        title='Right Hemisphere'
    )
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_clim(vmin, vmax)
    cbar = fig.colorbar(sm, ax=[ax1, ax2], shrink=0.6, location='right')
    cbar.set_label('Coupling Strength', fontsize=11)
    
    plt.suptitle(title, fontsize=14, y=0.95)
    plt.show()


# Main function for plotting Coupling Brain Map
def plot_all_unique_pairs(map_names, prepared_maps, name_dict):
    """
    Iterates through all unique pairs of maps, calculates Spearman Rho,
    determines direction (pos/neg), generates coupling map, and plots it.
    """
    unique_pairs = list(combinations(range(len(map_names)), 2))

    for i, (idx_a, idx_b) in enumerate(unique_pairs):
        name_a_raw = map_names[idx_a]
        name_b_raw = map_names[idx_b]
        data_a = prepared_maps[idx_a]
        data_b = prepared_maps[idx_b]
        
        # Calculate Spearman correlation for title & direction
        mask = (~np.isnan(data_a)) & (~np.isnan(data_b)) & (data_a != 0) & (data_b != 0)
        rho, _ = spearmanr(data_a[mask], data_b[mask])
        
        # If correlation is notably negative, flip the ranking logic
        # Threshold of -0.1 ensures weak negatives are treated as inverse
        is_positive = rho >= -0.1 
        
        # Generate the Coupling Map
        coupling_brain = generate_coupling_map(data_a, data_b, positive=is_positive)
        
        # Create Title
        title_a = get_readable_title(name_a_raw, name_dict)
        title_b = get_readable_title(name_b_raw, name_dict)
        
        full_title = f"{title_a} vs {title_b} (ρ={rho:.2f})"
    
        plot_coupling_map(coupling_brain, title=full_title)