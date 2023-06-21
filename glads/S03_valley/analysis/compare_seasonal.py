"""

Compare steady state mountain glacier simulations

"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib import gridspec
import netCDF4 as nc

import cmocean

import GladsPlot as gplt
from palettes.code import palettes

colors = np.array([[0.579, 0.677, 0.781, 1],
                   [0.199, 0.328, 0.492, 1],
                   [0.250, 0.250, 0.250, 1],
                   [0.929, 0.835, 0.408, 1],
                   [0.836, 0.590, 0.160, 1]])

rhow = 1000
g = 9.8
nu = 1.79e-6

def compare_seasonal(fnames, figname, labels, vmin=-1e6, vmax=1e6,
    Qmin=1e-3, Qmax=10, tslice=365+190+1, colors=colors):
    """Compare steady simulations"""

    n_cases = len(fnames)
    fig = plt.figure(figsize=(8, (1.5*n_cases)))
    hratios = 100*np.ones(n_cases+1)
    hratios[0] = 8
    gs = gridspec.GridSpec(n_cases+1, 4,
        width_ratios=[100, 2, 20, 100],
        height_ratios=hratios,
        left=0.1, right=0.9,
        top=0.9, bottom=0.1,
        wspace=0.1)
    axs = [fig.add_subplot(gs[i+1, 0]) for i in range(n_cases)]
    tax = fig.add_subplot(gs[1:, 3])
    cax_top = fig.add_subplot(gs[0, 0])
    cax_right = fig.add_subplot(gs[1:, 1])

    alphabet = ['a', 'b', 'c', 'd', 'e', 'f']
    for ii,fname in enumerate(fnames):
        with nc.Dataset(fname, 'r') as out:

            connect = out['connect'][:].data.T - 1
            connect_edge = out['connect_edge'][:].data.T.astype(int) - 1
            nodes = out['nodes'][:].data.T

            bed = np.vstack(out['bed'][:].data.T)
            thick = np.vstack(out['thick'][:].data.T)
            phi_z = rhow*g*bed

            h_sheet = out['h_sheet'][:].data.T
            phi = out['phi'][:].data.T
            N = out['N'][:].data.T

            pw = phi - phi_z
            pi = 910*g*thick
            ff = pw/pi

            Q = np.abs(out['Q'][:].data.T)
            S = out['S_channel'][:].data.T

            qxy = out['qs'][:].data.T
            qs = np.sqrt(qxy[:, 0]**2 + qxy[:, 1]**2)
            elements = out['elements'][:].data.T

            time = out['time'][:].data.T
            tt = time/86400/365
        
        mtri = Triangulation(nodes[:, 0], nodes[:, 1], connect)

        ax = axs[ii]
        pc = ax.tripcolor(mtri, pw[:, tslice], vmin=-vmax, vmax=vmax,
            cmap=cmocean.cm.curl)
        

        line_cmap = palettes.get_cmap('BrownYellow')
        lc = gplt.plot_edge_data(nodes, connect_edge, Q[:, tslice],
                    line_cmap, vmin=Qmin, vmax=Qmax)
        ax.add_collection(lc)

        ax.text(0.025, 1, alphabet[ii], fontweight='bold',
            transform=ax.transAxes)
        
        ax.text(0.95, 1, labels[ii], transform=ax.transAxes,
            ha='right')

        ax.set_aspect('equal')
        ax.set_xlim([0, 6e3])
        ax.set_ylim([-600, 600])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        tax.plot(tt-101, np.mean(ff, axis=0), color=colors[ii])
        tax.axvline((tt-101)[tslice], color='k')

    fig2, ax2 = plt.subplots()
    ax2.plot(tt-101, np.max(h_sheet, axis=0))

    cbar_top = fig.colorbar(pc, cax=cax_top, orientation='horizontal')
    cbar_right = fig.colorbar(lc, cax=cax_right)

    cax_top.xaxis.tick_top()
    cax_top.xaxis.set_label_position('top')

    cbar_top.set_label(r'$p_{\rm{w}}$ (Pa)')
    cbar_right.set_label(r'$Q~({\rm{m}}^3~{\rm{s}}^{-1})$')

    tax.set_xlim([0.25, 10/12])

    for ax in axs[:-1]:
        ax.set_xticklabels([])

    fig.savefig(figname)

if __name__=='__main__':
    cases = [101, 103, 105]
    labels = ['Turbulent 5/4', 'Laminar', 'Transition 3/2']
    fpattern = '../RUN/output_%03d_seasonal.nc'
    fnames = [fpattern % caseid for caseid in cases]
    figname = 'seasonal_valley_100.png'
    compare_seasonal(fnames, figname, labels,
        colors=colors[np.array([0, 2, 4]), :])

plt.show()
