"""

Plot 2x2 grid of open channel volume

    [Flat bed, Trough bed] x [SHMIP forcing, KAN_L forcing]

"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import netCDF4 as nc

import helpers
import defaults


figsize=(7, 4)


def plot_channel_volume(fnames, figname,
    x_band=defaults.x_bands[1], band_width=defaults.band_width, 
    figsize=figsize, labels=defaults.labels, 
    colors=defaults.colors, tlim=[1, 2],
    t_ticks=[1.0, 1.25, 1.5, 1.75, 2], ylim=None,
    t_ticklabels=None, xlabel=r'S (m$^2$)',
    rhow=1000, g=9.81, tslice=365+120):
    """
    Plot 2D floatation fraction maps and timeseries.

    [More description]

    Inputs:
    --------
    fnames : Iterable of str
        List of paths to model outputs for a simulation case
    
    figname : str
        Path to save figure
    
    Options:
    --------
    x_band : floats
        x distance in km for timeseries
    
    band_width : float
        Width of band for band-averaging
    
    figsize : (width, height) = (6, 7)
        Tuple in inches of figure size. 
    
    labels : list of str
        List of strings specifying labels for legend
    
    colors : (M, 3) or (M, 4) array-like of rgb or rgba values
    
    tlim : (tmin, tmax)
        x-axis bounds for timeseries panels
    
    t_ticks : Array-like
        x-axis ticks for timeseries panels

    ylim : (ymin, ymax) for ff panels or None
    
    Returns:
    --------
        fig : matplotlib figure object
    """

    # Start figure
    fig = plt.figure(figsize=figsize)

    gs = gridspec.GridSpec(2, 2, wspace=0.15, hspace=0.1,
        left=0.1, right=0.97, bottom=0.12, top=0.85)
    axs = np.array([[fig.add_subplot(gs[i,j]) for j in range(2)] for i in range(2)])

    dmesh = nc.Dataset('../glads/data/mesh/mesh_04.nc')

    for ii in range(len(fnames)):
        axii = axs.flat[ii]
        for jj in range(len(fnames[ii])):
            fname = fnames[ii][jj]
            print(fname)
            out = nc.Dataset(fname, 'r')
            nodes = out['nodes'][:].data.T

            # Get channel cross section
            S = out['S_channel'][:].data.T
            channel_volume = S*np.vstack(dmesh['tri/edge_length'][:].data)
            volume_slice = channel_volume[:, tslice]

            # axii.fill_between(tt, ff_lower, ff_upper, facecolor=colors[jj], alpha=0.3)
            # axii.plot(tt, sum_volume, color=colors[jj], label=labels[jj], linewidth=1)
            # axii.set_yscale('log')
            bins = np.linspace(-6, 3, 25)
            axii.hist(volume_slice, cumulative=False, log=True, facecolor=colors[jj], alpha=0.3,
                range=(10**bins[0], 10**bins[-1]), bins=10**bins, edgecolor=colors[jj],
                histtype='step', density=True)
            axii.set_xscale('log')

    axs.flat[0].legend(bbox_to_anchor=(0.2, 1.15, 1.8, 0.2), ncol=3, frameon=False)
    alphabet = ['a', 'b', 'c', 'd']
    for i,ax in enumerate(axs.flat):
        ax.grid()
        # ax.set_xlim(tlim)
        # ax.set_xticks(t_ticks)
        # 
        # if t_ticklabels:
        #     ax.set_xticklabels(t_ticklabels)
        # 
        # ax.text(0.05, 0.95, alphabet[i], va='top', transform=ax.transAxes,
        #     fontweight='bold')


    for ax in axs[0]:
        ax.set_xticklabels([])
    
    for ax in axs[:, 1]:
        ax.set_yticklabels([])

    for ax in axs[:, 0]:
        ax.set_ylabel(r'Count')
    
    for ax in axs[1]:
        ax.set_xlabel(xlabel)

    fig.savefig(figname, dpi=600)
    return fig

if __name__=='__main__':
    ## Global

    t_lim = [1 + 4/12, 1 + 10/12]
    t_ticks = 1 + np.arange(4, 11)/12
    t_ticklabels = ['4', '', '6', '', '8', '', '10']
    ff_ylim = [0, 1.75]

    cases = [[201, 202, 203, 204, 205],
         [101, 102, 103, 104, 105],
         [301, 302, 303, 304, 305],
         [101, 102, 103, 104, 105]
            ]

    patterns = ['../glads/00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc',
            '../glads/01_kan_l_forcing_shmip_topo/RUN/output_%03d_seasonal.nc',
            '../glads/02_shmip_forcing_synth_topo/RUN/output_%03d_seasonal.nc',
            '../glads/03_kan_l_forcing_synth_topo/RUN/output_%03d_seasonal.nc',
        ]

    fnames = [[patterns[j] % cases[j][i] for i in range(5)] for j in range(4)]
    figname = 'pressure_grid_new.png'
    fig_00 = plot_channel_volume(fnames, figname)

    plt.show()

