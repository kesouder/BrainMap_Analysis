import numpy as np
import matplotlib.pyplot as plt
from neuromaps import datasets, transforms, plotting, images

def plot_brain_map(map: dict, map_names: dict, brain_map_settings: dict):
    """
    Plots the brain map given in the map dictionary
    Parameters:
    map: the map as a dictionary with the needed parameters for feth_annotation
    map_names: dictionary of formal names for readability and plot titles
    Outputs: 
        Brain Maps!
    """
    map_paper, map_desc, map_space, map_den = map.values()
    #fetch source map and target map files
    src_map = datasets.fetch_annotation(**map)

    if map_desc == 'way100635':
        src_map = transforms.mni152_to_fslr(src_map, '32k')
        map_space = 'fsLR'
        map_den = '32k'
        
    settings = brain_map_settings.get(map_desc, {})
    cmap = settings.get('cmap', 'inferno')

    fig = plt.figure(figsize=(10, 4))
    fslr = datasets.fetch_atlas(map_space, map_den)
    surf_mesh_left = fslr['inflated'].L
    surf_mesh_right = fslr['inflated'].R
    data_full = images.load_data(src_map)

    if settings.get('vmin') ==  'special_perc' and settings.get('vmax') == 'special_perc':
        vmin, vmax = np.percentile(data_full[~np.isnan(data_full)], [10, 95])
    elif settings.get('vmin') is not None and settings.get('vmax') is not None:
        vmin, vmax = settings['vmin'], settings['vmax']
    else:
        vmin, vmax = np.percentile(data_full[~np.isnan(data_full)], [2, 98])

    data_l = images.load_data(src_map[0])
    data_r = images.load_data(src_map[1])
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    plotting.plot_surf(
        surf_mesh=surf_mesh_left,
        surf_map=data_l,
        hemi='left',
        view='lateral',
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        axes=ax1,
        title='Left hemisphere'
    )
    # right hemi
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    plotting.plot_surf(
        surf_mesh=surf_mesh_right,
        surf_map=data_r,
        hemi='right',
        view='lateral',
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        axes=ax2,
        title='Right hemisphere'
    )
    # color bar
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_clim(vmin, vmax)
    cbar = fig.colorbar(sm, ax=[ax1, ax2], shrink=0.6, location='right')
    cbar.set_label(f"{map_names.get(map_desc)}({map_space} {map_den})", fontsize=11)
    plt.suptitle(f"{map_names.get(map_desc)}", fontsize=14)
    plt.show()

    return None