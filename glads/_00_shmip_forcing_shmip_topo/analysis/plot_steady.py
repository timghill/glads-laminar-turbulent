"""
Plot floatation fraction for turbulent, laminar, and transition models

"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.tri import Triangulation
import netCDF4 as nc

import cmocean
from palettes.code import palettes, tools

import GladsPlot as gplt

# Define fnames
fname_pattern = '../RUN/output_%03d_steady.nc'
cases = [1, 2, 2]
n_cases = 3
figname = 'steady.png'

tslices = [-1]
n_times = len(tslices)

# Define the figure
fig = plt.figure(figsize=(8, 6))
gs = GridSpec(n_cases+2, 3,
    wspace=0.1, hspace=0.05,
    left=0.1, right=0.9, top=0.9, bottom=0.1,
    height_ratios=([8] + n_cases*[100] + [150]),
    width_ratios=(100, 100, 5))

axs = np.array([[fig.add_subplot(gs[i+1, j]) for j in range(2)] for i in range(n_cases)])
print(axs.shape)
alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']

axs_scatter = np.array([fig.add_subplot(gs[n_cases+1, j]) for j in range(2)])
labels = ['Turbulent', 'Laminar', 'Transition', 'Turbulent (3/2)', 'Transition (3/2)']

colors = np.array([[0.420, 0.510, 0.620, 1],
#                   [0.579, 0.677, 0.781, 1],
                   [0.500, 0.500, 0.500, 1],
                   [0.859, 0.683, 0.275, 1],
                    ])
#                   [0.929, 0.835, 0.408, 1]])

cm1 = palettes.get_cmap('BrownGray')
cm2 = palettes.get_cmap('blue-8').reversed()
cm = tools.join_cmaps(cm1, cm2, N1=50, N2=150, average=10)

# Start reading the data
for ii in range(n_cases):
    fname = fname_pattern % cases[ii]
    print(fname)

    out = nc.Dataset(fname, 'r')
    nodes = out['nodes'][:].data.T
    connect = out['connect'][:].data.T.astype(int) - 1
    connect_edge = out['connect_edge'][:].data.T.astype(int) - 1
    elements = out['elements'][:].data.T

    Q = np.abs(out['Q'][tslices, :].data.T)

    phi = out['phi'][tslices, :].data.T
    N = out['N'][tslices, :].data.T

    q = out['qs'][tslices, :].data.T
    qs = np.sqrt(q[:, 0]**2 + q[:, 1]**2)
    Re = qs/float(out['para/nu'][:].data)

    phi_0 = 9.81*1000*np.vstack(out['bed'][:].data)
    pw = phi - phi_0
    ff = pw/(N + pw)

    mtri = Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)

    ax = axs[ii, 0]
    fcolor = ax.tripcolor(mtri, ff[:, -1], cmap=cmocean.cm.dense, vmin=0, vmax=1)
    ax.set_aspect('equal')
    ax.set_xlim([0, 100])
    ax.set_ylim([0, 25])
    ax.set_yticks([0, 12.5, 25])

    lc = gplt.plot_edge_data(nodes/1e3, connect_edge, Q[:, -1],
        palettes.get_cmap('BrownYellow'), vmin=10, vmax=100)
    ax.add_collection(lc)
    ax.set_xticklabels([])
    ax.text(-0.05, 1.05, alphabet[ii*n_times], transform=ax.transAxes, fontweight='bold')

    ax = axs[ii, 1]
    rcolor = ax.tripcolor(mtri, Re[:, -1], cmap=cm, vmin=0, vmax=4000)
    ax.set_aspect('equal')
    ax.set_xlim([0, 100])
    ax.set_ylim([0, 25])
    ax.set_yticks([0, 12.5, 25])

    lc = gplt.plot_edge_data(nodes/1e3, connect_edge, Q[:, -1],
        palettes.get_cmap('BrownYellow'), vmin=10, vmax=100)
    ax.add_collection(lc)

    ax.set_xticklabels([])
    ax.set_yticklabels([])

        
#    if ii==0:
#        ax.text(95, 23, 't = %s a' % tlabels[jj], ha='right', color='w',
#                fontweight='bold', va='top')
#    if jj==1:
#        ax.text(95, 2, labels[ii], ha='right', color='w', fontweight='bold')

#    if jj==1:
#        ax.set_yticklabels([])

#    if ii<n_cases:
#        ax.set_xticklabels([])
    ax.text(-0.05, 1.05, alphabet[ii*n_times+1], transform=ax.transAxes, fontweight='bold')
    
    axs_scatter[0].scatter(nodes[:, 0]/1e3, ff[:, -1], 5, alpha=0.5, color=colors[ii],
            label=labels[ii])
    axs_scatter[0].set_ylim([0, 1])
    axs_scatter[0].set_xlim([0, 100])
    axs_scatter[0].grid()

    axs_scatter[1].scatter(elements[:, 0]/1e3, Re[:, -1], 5, alpha=0.5, color=colors[ii],
        label=labels[ii])
    axs_scatter[1].set_ylim([0, 4000])
    axs_scatter[1].set_xlim([0, 100])
    axs_scatter[1].grid()

axs_scatter[0].legend(markerscale=3)
axs_scatter[0].set_ylabel(r'$p_{\rm{w}}/p_{\rm{i}}$')

axs_scatter[0].text(-0.05, 1.05, 'e', transform=axs_scatter[0].transAxes, fontweight='bold')
axs_scatter[1].text(-0.05, 1.05, 'f', transform=axs_scatter[1].transAxes, fontweight='bold')

axs_scatter[0].set_xlabel('x (km)')
axs_scatter[1].set_xlabel('x (km)')

axs_scatter[1].set_ylabel('Re')

axs_scatter[1].yaxis.tick_right()
axs_scatter[1].yaxis.set_label_position('right')

cax1 = fig.add_subplot(gs[0, 0])
cax2 = fig.add_subplot(gs[0, 1])

fig.colorbar(fcolor, cax=cax1, orientation='horizontal')
fig.colorbar(rcolor, cax=cax2, orientation='horizontal', extend='max')

cax1.xaxis.tick_top()
cax1.xaxis.set_label_position('top')

cax2.xaxis.tick_top()
cax2.xaxis.set_label_position('top')

cax1.set_xlabel(r'$p_{\rm{w}}/p_{\rm{i}}$')
cax2.set_xlabel(r'$\rm{Re}$')

cax3 = fig.add_subplot(gs[1:4, 2])
cax3.set_yticks([10, 20, 40, 60, 80, 100])

cb3 = fig.colorbar(lc, cax=cax3, extend='both')
cb3.set_label(r'$Q~(\rm{m}^3~\rm{s}^{-1})$')


fig.savefig(figname, dpi=600)
plt.show()
