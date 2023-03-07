"""

Plot Reynolds number averaged within discrete bands

"""

import numpy as np
import netCDF4 as nc

from matplotlib import pyplot as plt
from matplotlib import gridspec


## CONFIG

# Define fnames
fname_pattern = '../RUN/output_%03d_seasonal.nc'
cases = [1, 2, 3, 4, 5]
n_cases = len(cases)
figname = 'floatation_composite.png'

# X bands for timeseries
x_bands = [15, 30, 70]
band_width = 3

# customizations
xlims = [1, 2]
xticks = [1, 1.25, 1.5, 1.75, 2]

## Start the figure
fig = plt.figure(figsize=(7, 6))

# Simple gridspec
gs = gridspec.GridSpec(len(x_bands), 2,
    left=0.1, right=0.95, wspace=0.2, hspace=0.1)
axs = np.array([[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(len(x_bands))])

# labels
alphabet = ['a', 'b', 'c', 'd', 'e', 'f']
labels = ['Turbulent 5/4', 'Turbulent 3/2', 'Laminar', 'Transition 5/4', 'Transition 3/2']

colors = np.array([[4.203e-01, 5.10e-01, 6.20e-01,1],
                   [0.5793, 0.677, 0.7812, 1],
                   [0.5, 0.5, 0.5, 1],
                   [0.859, 0.683, 0.275, 1],
                   [0.929, 0.835, 0.408, 1]])

rhow = 1000
g = 9.81
omega = 1/2000

# Make the plot
for ii in range(n_cases):
    fname = fname_pattern % cases[ii]
    print(fname)
    out = nc.Dataset(fname, 'r')

    time = out['time'][:].data
    tt = (time - time[0])/86400/365

    qxy = out['qs'][:].data.T
    nu = float(out['para/nu'][:])
    qs = np.sqrt(qxy[:, 0]**2 + qxy[:, 1]**2)
    Re = qs/nu

    phi = out['phi'][:].data.T
    N = out['N'][:].data.T

    bed = np.vstack(out['bed'][:].data)
    phi_elevation = rhow*g*bed

    pw = phi - phi_elevation
    ff = pw/(N + pw)

    nodes = out['nodes'][:].data.T
    elements = out['elements'][:].data.T

    # Timeseries
    for j, xb in enumerate(x_bands):
        axi = axs[j, 0]
        xmin = xb - band_width/2
        xmax = xb + band_width/2
        node_mask = np.logical_and(nodes[:, 0]/1e3>=xmin, nodes[:, 0]/1e3<xmax)
        f_mean = np.mean(ff[node_mask, :], axis=0)

        elem_mask = np.logical_and(elements[:, 0]/1e3>=xmin, elements[:, 0]/1e3<xmax)
        Re_mean = np.mean(Re[elem_mask, :], axis=0)

        axi.plot(tt, f_mean, label=labels[ii], color=colors[ii])
        axi.text(0.05, 0.95, alphabet[2*j], transform=axi.transAxes,
            fontweight='bold', va='top')

        axi.grid()
        axi.set_ylim([0, 1.2])
        axi.set_xlim(xlims)
        axi.set_xticks(xticks)

        axj = axs[j, 1]
        axj.plot(tt, omega*Re_mean, label=labels[ii], color=colors[ii])
        axj.text(0.05, 0.95, alphabet[2*j+1], transform=axj.transAxes,
            fontweight='bold', va='top')
        axj.grid()
        axj.set_xlim(xlims)
        axj.set_xticks(xticks)

    out.close()

axs[1, 0].set_ylabel(r'$p_{\rm{w}}/p_{\rm{i}}$')

axs[1, 1].set_ylabel(r'$\omega {\rm{Re}}$')

axs[0, 0].legend(ncol=3, bbox_to_anchor=(0.1, 1.05, 1.8, 0.2),
    frameon=False)

for ax in axs[:-1].flatten():
    ax.set_xticklabels([])

for j, xb in enumerate(x_bands):
    axi = axs[j, 1]
    axi.text(0.925, 0.925, 'x = %d km' % xb,
        transform=axi.transAxes, ha='right', va='top',
        backgroundcolor='w')

axs[-1, 0].set_xlabel('Year')
axs[-1, 1].set_xlabel('Year')

plt.show()


