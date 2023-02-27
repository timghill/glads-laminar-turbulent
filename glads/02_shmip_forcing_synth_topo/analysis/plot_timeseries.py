"""
Plot domain-averaged floatation fraction, space-time plots
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import netCDF4 as nc

import cmocean

fname_pattern = '../RUN/output_%03d_seasonal.nc'
cases = [911, 912, 913]
n_cases = len(cases)
figname = 'timeseries_91x.png'

# Make figure
fig = plt.figure(figsize=(8, 3.5))
gs = GridSpec(1, n_cases+1, width_ratios=(100, 100, 100, 8),
    left=0.075, right=0.925, bottom=0.175)

axs = np.array([fig.add_subplot(gs[0, ii]) for ii in range(n_cases)])

# Define grid for average
xedge = np.arange(0, 101, 1)
xmid = 0.5*(xedge[:-1] + xedge[1:])

# Read and plot each case
for jj in range(n_cases):
    fname = fname_pattern % cases[jj]
    out = nc.Dataset(fname)

    pw = out['phi'][:].data.T
    N = out['N'][:].data.T
    ff = pw/(N + pw)

    nodes = out['nodes'][:].data.T/1e3
    
    times = out['time'][:].data/86400/365 - 100

    [xx, yy] = np.meshgrid(xmid, times)
    print(xx.shape)
    ff_grid = np.zeros((len(xmid), ff.shape[1]))

    for kk in range(len(xmid)):
        x_mask = np.logical_and(nodes[:, 0]>=xedge[kk], nodes[:, 0]<xedge[kk+1])
        ff_grid[kk, :] = np.mean(ff[x_mask, :], axis=0)
    
    print(ff_grid.shape)
    ax = axs[jj]
    pcolor = ax.pcolormesh(xx, yy, ff_grid.T, cmap=cmocean.cm.dense, vmin=0, vmax=1)
    ax.contour(xx, yy, ff_grid.T, levels=[1], colors='w', linewidths=0.5)
    if jj>0:
        ax.set_yticklabels([])

axs[0].set_ylabel('Year')
axs[1].set_xlabel('x (km)')
cax = fig.add_subplot(gs[0, -1])
cb = fig.colorbar(pcolor, cax=cax)
cb.set_label(r'$p_{\rm{w}}/p_{\rm{i}}$')

fig.savefig(figname, dpi=600)
plt.show()
