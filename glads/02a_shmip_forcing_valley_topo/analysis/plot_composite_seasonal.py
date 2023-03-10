"""

Plot timeseries for discrete elevation bands and 2D snapshot maps

"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.tri import Triangulation
import netCDF4 as nc

import cmocean
from palettes.code import palettes, tools

import GladsPlot as gplt

## CONFIG

# Define fnames
fname_pattern = '../RUN/output_%03d_seasonal.nc'
cases = [301, 302, 303, 304, 305]
n_cases = len(cases)
figname = 'floatation_composite.png'

# Time slice for 2D snapshots
tslice = 180 + 365

# X bands for timeseries
x_bands = [15, 30, 70]
band_width = 3

## Start the figure
fig = plt.figure(figsize=(7, 6))

# A global gridspec giving two columns to work with
global_gs = gridspec.GridSpec(1, 3, wspace=0.05, hspace=0.2, 
    width_ratios = (100, 30, 115), left=0.1, right=0.98, bottom=0.08)

# Left column: 3 timeseries panels
gs_timeseries = gridspec.GridSpecFromSubplotSpec(len(x_bands), 1, 
    subplot_spec=global_gs[:, 2])

# Right column: 5 maps with space for colorbars
hratios = 100*np.ones(n_cases+2)
hratios[0] = 8
hratios[-1] = 150
gs_maps = gridspec.GridSpecFromSubplotSpec(n_cases+2, 2, 
    subplot_spec=global_gs[:, 0], width_ratios=(100, 5), height_ratios=hratios,
    hspace=0.04, wspace=0.1)

print(gs_maps)

# Initialize axes
axs_timeseries = np.array([fig.add_subplot(gs_timeseries[i, 0]) for i in range(len(x_bands))])
axs_maps = np.array([fig.add_subplot(gs_maps[i+1, 0]) for i in range(n_cases)])
ax_scatter = fig.add_subplot(gs_maps[-1, 0])

labels = ['Turbulent 5/4', 'Turbulent 3/2', 'Laminar', 'Transition 5/4', 'Transition 3/2']

# Set style for panel labels
time_alphabet = ['g', 'h', 'i']
map_alphabet = ['a', 'b', 'c', 'd', 'e', 'f']
text_args = {'fontweight':'bold'}

# Set colours for 5 models
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
print(colors)

colors = np.array([[4.203e-01, 5.10e-01, 6.20e-01,1],
                   [0.5793, 0.677, 0.7812, 1],
                   [0.5, 0.5, 0.5, 1],
                   [0.859, 0.683, 0.275, 1],
                   [0.929, 0.835, 0.408, 1]])

colors = colors[:, :3]
print(colors)

# Start reading the data
for ii in range(n_cases):
    fname = fname_pattern % cases[ii]
    print(fname)

    out = nc.Dataset(fname, 'r')
    nodes = out['nodes'][:].data.T
    connect = out['connect'][:].data.T.astype(int) - 1
    connect_edge = out['connect_edge'][:].data.T.astype(int) - 1

    # Channel fields
    Q = np.abs(out['Q'][:, :].data.T)

    # Get floatation fraction
    phi = out['phi'][:, :].data.T
    N = out['N'][:, :].data.T
    phi_0 = 9.81*1000*np.vstack(out['bed'][:].data)
    pw = phi - phi_0
    ff = pw/(N + pw)

    tt = out['time'][:].data/86400/365 - 100

    # Initialize triangulation for faster plotting
    mtri = Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)
    
    # Map panel
    mapax = axs_maps[ii]
    fcolor = mapax.tripcolor(mtri, ff[:, tslice], cmap=cmocean.cm.dense, vmin=0, vmax=1)
    mapax.set_aspect('equal')
    mapax.set_xlim([0, 100])
    mapax.set_ylim([0, 25])
    mapax.set_yticks([0, 12.5, 25])

    lc = gplt.plot_edge_data(nodes/1e3, connect_edge, Q[:, tslice],
        palettes.get_cmap('BrownYellow'), vmin=10, vmax=100)
    mapax.add_collection(lc)
    
    # if ii>0:
        # mapax.set_yticklabels([])

    # mapax.yaxis.tick_right()
    # mapax.yaxis.set_label_position('right')

    if ii<n_cases:
        mapax.set_xticklabels([])
    
    mapax.text(-0.25, 0.95, map_alphabet[ii], transform=mapax.transAxes, **text_args)
    print(ff.shape)
    ax_scatter.scatter(nodes[:, 0]/1e3, ff[:, tslice], 5, alpha=0.5, color=colors[ii],
        label=labels[ii])
    ax_scatter.set_ylim([0, 1.25])
    ax_scatter.set_xlim([0, 100])
    ax_scatter.grid()

    if ii==0:
        ax_scatter.set_xlabel('x (km)')

    # Timeseries
    for j, xb in enumerate(x_bands):
        xmin = xb - band_width/2
        xmax = xb + band_width/2
        node_mask = np.logical_and(nodes[:, 0]/1e3>=xmin, nodes[:, 0]/1e3<xmax)
        f_mean = np.mean(ff[node_mask, :], axis=0)
        timeax = axs_timeseries[j]
        timeax.plot(tt, f_mean, label=labels[ii], color=colors[ii])

        mapax.axvline(xb, color='w', linewidth=0.5)
        timeax.axvline(tslice/365, color='k', linewidth=0.5)

        timeax.text(0.05, 0.9, time_alphabet[j], transform=timeax.transAxes, **text_args)


# ax_scatter.legend(markerscale=3)
ax_scatter.set_ylabel(r'$p_{\rm{w}}/p_{\rm{i}}$')
# ax_scatter.yaxis.tick_right()
# ax_scatter.yaxis.set_label_position('right')
ax_scatter.text(0.05, 0.85, map_alphabet[n_cases], transform=ax_scatter.transAxes, **text_args)

axs_maps[2].set_ylabel('y (km)')

cax1 = fig.add_subplot(gs_maps[1:6, 1])
cax2 = fig.add_subplot(gs_maps[0, 0])

fig.colorbar(fcolor, cax=cax1, extend='max')
fig.colorbar(lc, cax=cax2, orientation='horizontal', extend='both')

cax1.xaxis.tick_top()
cax1.xaxis.set_label_position('top')

cax2.set_xticks([10, 20, 40, 60, 80, 100])

cax2.xaxis.tick_top()
cax2.xaxis.set_label_position('top')

# cax1.set_ylabel(r'$p_{\rm{w}}/p_{\rm{i}}$')
cax1.text(0, 1.1, r'$p_{\rm{w}}/p_{\rm{i}}$', transform=cax1.transAxes)
cax2.set_xlabel(r'$Q~(\rm{m}^3~\rm{s}^{-1})$')

axs_timeseries[0].legend(bbox_to_anchor=[0, 1.02, 1., 0.102], loc='lower left',
    ncol=2, mode='expand', borderaxespad=0.05, frameon=False, borderpad=0)
for j, xb in enumerate(x_bands):
    axi = axs_timeseries[j]
    if j<len(x_bands)-1:
        axi.set_xticklabels([])
    axi.set_xlim([1, 2])
    axi.set_ylim([0, 1.5])
    axi.set_xticks([1.0, 1.25, 1.5, 1.75, 2])
    axi.grid()

axs_timeseries[-1].set_xlabel('Time (a)')
axs_timeseries[1].set_ylabel(r'$p_{\rm{w}}/p_{\rm{i}}$')

fig.savefig(figname, dpi=600)
plt.show()
