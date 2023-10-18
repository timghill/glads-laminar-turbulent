"""

Plot floatation fraction maps and timeseries.

Function plot_pressure_maps_timeseries does the work,
this can be called from external scripts

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
import defaults

zorders = np.array([2, 2, 3, 2, 2]) + 2


figsize=(5, 5)

def plot_timeseries(fnames, figname, tslice=defaults.tslice, 
    x_bands=[30], band_width=defaults.band_width, 
    figsize=figsize,labels=defaults.labels, 
    colors=defaults.colors, map_cmap=defaults.cmaps['floatation'],
    line_cmap=defaults.cmaps['Q'], Qmin=10, Qmax=100,
    t_lim=[1, 2], t_ticks=[1.0, 1.25, 1.5, 1.75, 2], ff_ylim=[0, 1.5],
    t_ticklabels=None, t_xlabel='Year', ff_yticks=[0, 0.25, 0.5, 0.75, 1, 1.25, 1.5],
    melt_forcing='SHMIP', fill_between=False,
    lws=defaults.linewidths, linestyles=defaults.linestyles,
    zorders=zorders, prefix='../../',
    boxes=[(1+4/12, 1+5/12), (1+5/12, 1+6/12), (1+6/12, 1+8/12)],
    boxmonths=[('June', '', '', '', 'July'), ('May','','','', 'June'),
                ('July', '', 'Aug', '', 'Sep')]):
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
        for xb in x_bands.Temp 

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

    t_ticks : Array-like1
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
        tt_temp = np.loadtxt(prefix + 'glads/data/kan_l_melt/KAN_L_2014_temp_clipped.txt', delimiter=',')
        tt_days = tt_temp[:, 0]
        temp_sl = tt_temp[:, 1]
        lr = -0.005
        DT = lr*390
        temp_fun = lambda t: np.maximum(0*t, DT + np.interp(t%1, tt_days/365, temp_sl, left=0, right=0))
    elif melt_forcing=='KANadj':
        tt_temp = np.loadtxt(prefix + 'glads/data/kan_l_melt/KAN_L_2014_temp_adjusted.txt', delimiter=',')
        tt_days = tt_temp[:, 0]
        temp_sl = tt_temp[:, 1]
        lr = -0.0075
        DT = lr*390
        temp_fun = lambda t: np.maximum(0*t, DT + np.interp(t%1, tt_days/365, temp_sl, left=0, right=0))
    elif melt_forcing=='SHMIP':
        temp_fun = lambda t: np.maximum(-16*np.cos(2*np.pi*t) - 5, 0*t)
    elif melt_forcing=='SHMIPadj':
        a = 9.0684
        DT = 390*0.0075
        DT_term = 0.005*390
        temp_fun = lambda t: -a*np.cos(2*np.pi*t) + a*(DT - 5)/16 - DT_term


    ## Start the figure
    fig = plt.figure(figsize=figsize)
    figs_focus = [plt.figure(figsize=(4, 2.5)) for ii in range(len(boxes))]
    axs_focus = [ff.add_subplot() for ff in figs_focus]

    # A global gridspec giving two columns to work with
    # global_gs = gridspec.GridSpec(1, 3, **gs_kwargs)

    # Left column: 3 timeseries panels and melt forcing
    hratios = 100*np.ones(len(x_bands)+2)
    hratios[0] = 8
    hratios[2] = 67
    gs_timeseries = gridspec.GridSpec(len(x_bands) + 2, 1, 
        height_ratios=hratios, left=0.125, right=0.95,
        bottom=0.1, top=0.925)

    # Initialize axes
    axs_timeseries = np.array([fig.add_subplot(gs_timeseries[i+1, 0]) for i in range(len(x_bands) + 1)])

    # Set style for panel labels
    time_alphabet = ['a', 'b', 'c', 'd']
    text_args = {'fontweight':'bold'}
    # Start reading the data
    for ii in range(n_cases):
        fname = fnames[ii]
        print(fname)

        with nc.Dataset(fname, 'r') as out:
            walltime = float(out['model/wallclock'][:])
            # print('Walltime (hours):', walltime/3600)

            nodes = out['nodes'][:].data.T
            connect = out['connect'][:].data.T.astype(int) - 1
            connect_edge = out['connect_edge'][:].data.T.astype(int) - 1

            # Channel fields
            Q = np.abs(out['Q'][:, :].data.T)

            # Get floatation fraction.lege
            phi = out['phi'][:, :].data.T
            N = out['N'][:, :].data.T
            phi_0 = 9.81*1000*np.vstack(out['bed'][:].data)
            pw = phi - phi_0
            ff = pw/(N + pw)

            tt = out['time'][:].data/86400/365 - 100

        with nc.Dataset(prefix + 'glads/data/mesh/mesh_04.nc', 'r') as dmesh:
            node_area = dmesh['tri/area_nodes'][:].data

        # Timeseries
        for j, xb in enumerate(x_bands):
            xmin = xb - band_width/2
            xmax = xb + band_width/2
            node_mask = np.logical_and(nodes[:, 0]/1e3>=xmin, nodes[:, 0]/1e3<xmax)
            f_mean = np.sum(ff[node_mask, :]*np.vstack(node_area[node_mask]), axis=0)/np.sum(node_area[node_mask])
            f_lower = np.quantile(ff[node_mask, :], 0.025, axis=0)
            f_upper = np.quantile(ff[node_mask, :], 0.975, axis=0)
            timeax = axs_timeseries[j]

            timeax.plot(tt, f_mean, label=labels[ii],
                color=colors[ii], linewidth=lws[ii], zorder=zorders[ii],
                linestyle=linestyles[ii])

            # mapax.axvline(xb, color='w', linewidth=0.5)
            # timeax.axvline(tslice/365, color='k', linewidth=0.5)

            timeax.text(0.025, 0.95, time_alphabet[j], transform=timeax.transAxes,
            va='top', ha='left', **text_args)

            timeax.text(0.95, 0.95, '%d km' % xb, transform=timeax.transAxes,
            va='top', ha='right', **text_args)

            # timeax.axhline(1.0, color='grey', linewidth=0.5, linestyle='solid', zorder=2)

            timeax.set_xticklabels([])

            for j in range(len(boxes)):
                fbox = figs_focus[j]
                abox = axs_focus[j]
                abox.plot(tt, f_mean, label=labels[ii],
                    color=colors[ii], linewidth=lws[ii], zorder=zorders[ii],
                    linestyle=linestyles[ii])

    melt = temp_fun(tt)
    melt_ax = axs_timeseries[-1]
    melt_ax.plot(tt, melt, color='k', linewidth=1)

    melt_ax.set_ylabel('Temperature ($^\circ$C)')
    melt_ax.set_ylim([0, 12])
    melt_ax.set_yticks([0, 4, 8, 12])
    melt_ax.grid(linestyle=':', linewidth=0.5)
    # melt_ax.axvline(tslice/365, color='k', linewidth=0.5)
    melt_ax.text(0.025, 0.95, time_alphabet[j+1], transform=melt_ax.transAxes,
        va='top', ha='left', **text_args)

    axs_timeseries[0].legend(bbox_to_anchor=[0, 1.02, 1., 0.102], loc='lower left',
        ncol=2, mode='expand', borderaxespad=0.05, frameon=False, borderpad=0)
    for j, xb in enumerate(x_bands):
        axi = axs_timeseries[j]
        axi.set_xticklabels([])
        axi.set_xlim(t_lim)
        axi.set_ylim(ff_ylim)
        axi.set_yticks(ff_yticks)
        axi.set_xticks(t_ticks)
        axi.grid(linestyle=':', linewidth=0.5)
    

    melt_ax.set_xlim(t_lim)
    melt_ax.set_xticks(t_ticks)
    if t_ticklabels:
        melt_ax.set_xticklabels(t_ticklabels)

    axs_timeseries[-1].set_xlabel(t_xlabel)
    axs_timeseries[0].set_ylabel(r'$p_{\rm{w}}/p_{\rm{i}}$')
    
    fig.savefig(figname, dpi=300)
    for j in range(len(boxes)):
        bbox = boxes[j]
        axi = axs_focus[j]
        # axi.set_xticklabels([])
        axi.set_xlim(bbox)
        axi.set_xticks(np.linspace(bbox[0], bbox[1], 5))
        axi.set_xticklabels(boxmonths[j])
        axi.set_ylim([0.2, 1.6])
        axi.grid(linestyle=':', linewidth=0.5)

        axi.set_ylabel(r'$p_{\rm{w}}/p_{\rm{i}}$')
        # axi.set_facecolor((0.85, 0.85, 0.85))

        fc = axs_timeseries[0].fill_betweenx([0,1.5], [bbox[0], bbox[0]], [bbox[1], bbox[1]],
            color=(0.85, 0.85, 0.85), alpha=1, zorder=0)

        fc2 = melt_ax.fill_betweenx([0,12], [bbox[0], bbox[0]], [bbox[1], bbox[1]],
            color=(0.85, 0.85, 0.85), alpha=1, zorder=0)

        fig.savefig('timeseries_%02d.png' % j, dpi=400)
        fc.set_color((1, 1, 1, 0))
        fc2.set_color((1,1,1,0))

        figs_focus[j].subplots_adjust(left=0.15, right=0.95, top=0.95)
        figs_focus[j].savefig('timeseries_focus_%02d.png' % j, dpi=400)

    fig.savefig(figname, dpi=400)

    return fig

if __name__=='__main__':
    cases = [1, 2, 3, 4, 5]
    fpattern = '../../glads/01_kan_forcing/RUN/output_%03d_seasonal.nc'
    fnames = [fpattern % ii for ii in cases]
    figname = 'timeseries.png'
    fig = plot_timeseries(fnames, figname, tslice=365+190,
        melt_forcing='KAN', t_ticks=[1 + 4/12, 1 + 6/12, 1 + 8/12, 1+10/12],
        t_ticklabels=['May', 'July', 'Sep', 'Nov'], t_lim=[1+4/12, 1+10/12],
        t_xlabel='Month')
    plt.show()
