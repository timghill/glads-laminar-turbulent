"""
Compare ISSM and MATLAB model flotation fraction and channel area
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.tri import Triangulation
import netCDF4 as nc

import cmocean

import defaults
import GladsPlot as gplt

models = ['Turbulent 5/4', 'Turbulent 3/2', 'Laminar',
            'Transition 5/4', 'Transition 3/2']

def plot_issm(mat_fnames, issm_fnames, figname,
    colors=defaults.colors, xb=30e3, bw=5e3, Qb=15e3,
    models=models,
    tslice=220, map_cmap=defaults.cmaps['floatation'],
    line_cmap=defaults.cmaps['Q'],
    Qmin=1, Qmax=100,
    ):
    # Set up the figure and axes
    fig = plt.figure(figsize=(7, 6))
    gs = GridSpec(3, 2, left=0.09, right=0.875, bottom=0.125, top=0.95,
        wspace=0.1, hspace=0.1,
        width_ratios=(100, 5))
    axs = np.array([fig.add_subplot(gs[i,0]) for i in range(3)])
    cax1 = fig.add_subplot(gs[:2, 1])
    cax2 = fig.add_subplot(gs[-1, 1])
    alphabet = ['a', 'b', 'c', 'd', 'e']
    with nc.Dataset('../../../glads/data/mesh/mesh_04.nc', 'r') as dmesh:
        edge_length = np.vstack(dmesh['tri/edge_length'][:].data)
        nodes = dmesh['tri/nodes'][:].data.T
        connect = dmesh['tri/connect'][:].data.T.astype(int) - 1
    
    mtri = Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)

    ii=4
    mfname = mat_fnames[ii]
    ifname = issm_fnames[ii]

    with nc.Dataset(mfname, 'r') as mdata:
        nodes = mdata['nodes'][:].data.T
        phi = mdata['phi'][:].data.T
        N = mdata['N'][:].data.T
        phi_0 = 9.8*1000*np.vstack(mdata['bed'][:].data)
        pw = phi - phi_0
        ff = pw/(N + pw)
        mtt = mdata['time'][:].data/86400/365 - 100
        Q = np.abs(mdata['Q'][:, :].data.T)
        connect_edge = mdata['connect_edge'][:].data.T.astype(int) - 1
    
    with nc.Dataset(ifname, 'r') as idata:
        N_issm = idata['N'][:].data.T
        phi_issm = idata['phi'][:].data.T
        pw_issm = phi_issm - np.vstack(idata['phi_bed'][:].data)
        ff_issm = pw_issm/(N_issm + pw_issm)
        Q_issm = np.abs(idata['Q'][:, :].data.T)
        itt = idata['time'][:].data.T

    print('ISSM shape:')
    print(ff_issm.shape)

    print('MATLAB shape:')
    print(ff.shape)
    itt = itt - 8
    tslice_issm = int(9*(365/5) + tslice/5)
    tslice_mat = tslice + 365

    print('MATLAB time')
    print(mtt[tslice_mat])

    print('ISSM time')
    print(itt[tslice_issm])

    axm = axs[0]
    # ax.plot(mtt, np.mean(ff[nodemask, :], axis=0),
    #     color=colors[ii])
    # ax.plot(itt, np.mean(ff_issm[nodemask, :], axis=0),
    #     color=colors[ii], linestyle=':')
    # ax.text(0.05, 0.95, alphabet[ii], va='top', fontweight='bold', transform=ax.transAxes)
    # ax.text(0.95, 0.95, models[ii], va='top', ha='right', transform=ax.transAxes)
    # ax.grid(linewidth=0.5, linestyle=':')
    mat_pc = axm.tripcolor(mtri, ff[:, tslice_mat],
        cmap=map_cmap, vmin=0, vmax=1)
    lcm = gplt.plot_edge_data(nodes/1e3, connect_edge, Q[:, tslice_mat],
            line_cmap, vmin=Qmin, vmax=Qmax)
    axm.add_collection(lcm)
    
    axi = axs[1]
    issm_pc = axi.tripcolor(mtri, ff_issm[:, tslice_issm],
        cmap=map_cmap, vmin=0, vmax=1)
    
    lci = gplt.plot_edge_data(nodes/1e3, connect_edge, Q_issm[:, tslice_issm],
            line_cmap, vmin=Qmin, vmax=Qmax)
    axi.add_collection(lci)
    # axs.flat[-1].set_visible(False)

    axd = axs[2]
    diff_pc = axd.tripcolor(mtri, ff_issm[:, tslice_issm] - ff[:, tslice],
        cmap=cmocean.cm.balance, vmin=-0.1, vmax=0.1)

    mat_cbar = fig.colorbar(mat_pc, cax=cax1)
    mat_cbar.set_label('Flotation fraction')

    diff_cbar = fig.colorbar(diff_pc, cax=cax2)
    diff_cbar.set_label(r'$\Delta$ Flotation fraction')
    """
    for ax in axs.flat:
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    for ax in axs[0, :2]:
        ax.set_xticklabels([])
    
    for ax in axs[:, 1:].flat:
        ax.set_yticklabels([])
    
    

    fig.text(0.01, 0.5, r'$p_{\rm{w}}/p_{\rm{i}}$', va='center', rotation='vertical')
    fig.text(0.5, 0.015, 'Month', va='bottom')
    axs.flat[4].legend(labels=('Matlab', 'ISSM'),
        bbox_to_anchor=(1.1, 0, 1, 0.8), mode='expand',
        frameon=False)
    """

    for ax in axs:
        ax.set_aspect('equal')
        ax.set_xlim([0, 100])
        ax.set_ylim([0, 25])

    axs[0].set_title('Day %d' % tslice)
    if figname:
        fig.savefig(figname, dpi=400)
        

if __name__ == '__main__':
    prefix = '/home/tghill/projects/def-gflowers/tghill/laminar-turbulent/'
    mat_fnames = [prefix + 'glads/00_synth_forcing/RUN/output_%03d_seasonal.nc'%i for i in range(1,6)]
    issm_fnames = [prefix + 'issm/00_synth_forcing/RUN/output_%03d.nc'%i for i in range(1,6)]
    tslice = 160
    figname = 'glads_issm_spatial_%d.png' % tslice
    plot_issm(mat_fnames, issm_fnames, figname, tslice=tslice)

    plt.show()
