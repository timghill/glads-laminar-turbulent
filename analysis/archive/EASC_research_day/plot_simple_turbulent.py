"""

Simple timeseries to compare with Andrews et al.

"""

"""

Plot floatation fraction maps and timeseries

"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.tri import Triangulation
import netCDF4 as nc

import cmocean
from palettes.code import palettes, tools

import GladsPlot as gplt

import helpers
# import defaults

figsize=(5, 4)

gs_kwargs=dict(wspace=0.05, hspace=0.2, 
        height_ratios=(2, 1), top=0.9)

labels = ['Turbulent', 'Laminar', 'Transition']
colors = np.array([[0.420, 0.510, 0.620, 1],
                #    [0.579, 0.677, 0.781, 1],
                   [0.500, 0.500, 0.500, 1],
                   [0.859, 0.683, 0.275, 1],
                #    [0.929, 0.835, 0.408, 1]]
                ])

def plot_pressure_maps_timeseries(fnames, figname, tslice=182, 
    xb=20, band_width=5, 
    figsize=figsize, gs_kwargs=gs_kwargs, labels=labels, 
    colors=colors, map_cmap=cmocean.cm.dense,
    line_cmap=palettes.get_cmap('BrownYellow'), Qmin=10, Qmax=100,
    t_lim=[1, 2], t_ticks=[1.0, 1.25, 1.5, 1.75, 2], ff_ylim=[0, 1.5],
    t_ticklabels=None, t_xlabel='Year',
    melt_forcing='SHMIP', fill_between=False):
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
    tslice : int
        Time index for 2D pressure maps
    
    x_bands : Array-like of floats
        x distances in km for timeseries. Timeseries are computed for mean
        floatation fraction in bands [xb-band_width/2, xb+band_width/2]
        for xb in x_bands.
    
    band_width : float
        Width of bands for band-averaging
    
    figsize : (width, height) = (6, 7)
        Tuple in inches of figure size. 
    
    gs_Kwargs : dict
        Dictionary of options passed to gridspec.GridSpec for global config
    
    labels : list of str
        List of strings specifying labels for legend
    
    colors : (M, 3) or (M, 4) array-like of rgb or rgba values

    map_cmap : LinearSegmentedColormap 
        Colormap object for tripcolor panels
    
    line_cmap : LinearSegmentedColormap
        Colormap object for plotting channel discharge
    
    Qmin, Qmax : float
        Min and max discharge for channel discharge colorbar
    
    tlim : (tmin, tmax)
        x-axis bounds for timeseries panels
    
    t_ticks : Array-like
        x-axis ticks for timeseries panels

    ff_ylim : (ymin, ymax) for ff panels
    
    Returns:
    --------
        fig : matplotlib figure object
    
    """ 
    ## CONFIG

    n_cases = len(fnames)

    # Sort out melt forcing
    if melt_forcing=='KAN':
        tt_temp = np.loadtxt('../../glads/data/kan_l_melt/KAN_L_2014_temp_clipped.txt', delimiter=',')
        tt_days = tt_temp[:, 0]
        temp_sl = tt_temp[:, 1]
        lr = -0.005
        DT = lr*(350 + 740)
        temp_fun = lambda t: np.maximum(-1000*t, DT + np.interp(t%1, tt_days/365, temp_sl, left=0, right=0))
    elif melt_forcing=='SHMIP':
        temp_fun = lambda t: np.maximum(-16*np.cos(2*np.pi*t) - 5, 0*t)
    elif melt_forcing=='SHMIPadj':
        a = 9.0684
        DT = 390*0.0075
        temp_fun = lambda t: -a*np.cos(2*np.pi*t) + a*(DT - 5)/16


    ## Start the figure
    fig = plt.figure(figsize=figsize)

    gs = gridspec.GridSpec(2, 1, **gs_kwargs)

    axs = np.array([fig.add_subplot(gs[i, 0]) for i in range(2)])

    alphabet = ['a', 'b']
    text_args = {'fontweight':'bold'}

    # Start reading the data
    for ii in range(n_cases):
        fname = fnames[ii]
        print(fname)

        out = nc.Dataset(fname, 'r')

        nodes = out['nodes'][:].data.T
        connect = out['connect'][:].data.T.astype(int) - 1
        connect_edge = out['connect_edge'][:].data.T.astype(int) - 1

        # Get floatation fraction
        phi = out['phi'][:, :].data.T
        N = out['N'][:, :].data.T
        phi_0 = 9.81*1000*np.vstack(out['bed'][:].data)
        pw = phi - phi_0
        ff = pw/(N + pw)

        tt = out['time'][:].data/86400/365 - 100

        # for j, xb in enumerate(x_bands):
        xmin = xb - band_width/2
        xmax = xb + band_width/2
        node_mask = np.logical_and(nodes[:, 0]/1e3>=xmin, nodes[:, 0]/1e3<xmax)
        f_mean = np.mean(ff[node_mask, :], axis=0)
        f_lower = np.quantile(ff[node_mask, :], 0.025, axis=0)
        f_upper = np.quantile(ff[node_mask, :], 0.975, axis=0)
        timeax = axs[0]
        timeax.plot(tt*365, f_mean, label=labels[ii], color=colors[ii])#, linewidth=1)

        timeax.text(0.025, 0.95, alphabet[0], transform=timeax.transAxes,
            va='top', ha='left', **text_args)

        timeax.set_ylim([0, 2])

        out.close()

    melt = temp_fun(tt)
    ax_right = axs[-1]
    ax_right.plot(tt*365, melt, color='k')
    
    ax_right.set_ylabel(r'Temp ($^\circ{\rm{C}}$)')
    ax_right.set_ylim([-15, 10])
    ax_right.set_yticks([-15, -10, -5, 0, 5, 10])

    ax_right.text(0.025, 0.95, alphabet[1], transform=ax_right.transAxes,
            va='top', ha='left', **text_args)
    # axs[0].legend(bbox_to_anchor=[0.1, 1.04, 0.9, 0.102], loc='lower left',
    #     ncol=3, mode='expand', borderaxespad=0.05, frameon=False, borderpad=0)

    for j in range(2):
        axi = axs[j]

        axi.grid()
        axi.set_xlim(t_lim)

    axs[0].set_xticklabels([])
    axs[0].set_ylabel('Floatation fraction')
    axs[1].set_xlabel('Day of 2014')
    fig.savefig(figname, dpi=600)
    return fig


if __name__ == '__main__':
    # t_ticks = [1 + 4/12, 1 + 6/12, 1 + 8/12, 1 + 10/12]
    # t_ticklabels = ['4', '6', '8', '10']
    # t_lim = [t_ticks[0], t_ticks[-1]]
    t_lim = [140, 265]
    t_xlabel = 'Day of year'

    ## Case 01: Flat topo, KAN_L forcing
    KAN_tslice = 569
    # cases = [2, 3, 5]
    cases = [102]
    pattern = '../../glads/01_kan_l_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'
    # pattern = '../../glads/02a_shmip_forcing_valley_topo/RUN/output_%03d_seasonal.nc'
    fnames = [pattern % caseid for caseid in cases]
    figname = '01_pressure_turb.png'
    fig_01 = plot_pressure_maps_timeseries(fnames, figname, tslice=KAN_tslice, Qmin=1, melt_forcing='KAN',
        t_xlabel=t_xlabel, t_lim=t_lim)
    plt.show()
