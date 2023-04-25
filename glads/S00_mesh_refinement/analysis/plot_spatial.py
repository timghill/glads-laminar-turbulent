"""
Simple plot of floatation fraction and channel discharge
for steady simulations
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation

import netCDF4 as nc

import cmocean

import GladsPlot as gplt
from palettes.code import palettes, tools

def plot_steady(fname, figname):
    fig, ax = plt.subplots(figsize=(6, 2.5))

    out = nc.Dataset(fname, 'r')
    phi = out['phi'][-1].data
    rhow = 1000
    g = 9.8
    phi0 = rhow*g*out['bed'][:].data
    pw = phi - phi0

    N = out['N'][-1].data
    pi = N + pw
    ff = pw/pi

    Q = out['Q'][-1].data.T

    nodes = out['nodes'][:].data.T
    connect = out['connect'][:].data.T.astype(int) - 1
    connect_edge = out['connect_edge'][:].data.T.astype(int) - 1
    mtri = Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)
    
    # Plot floatation fraction
    fcolor = ax.tripcolor(mtri, ff, cmap=cmocean.cm.dense, vmin=0, vmax=1)

    # Color edges by discharge
    lc = gplt.plot_edge_data(nodes/1e3, connect_edge, np.abs(Q),
        palettes.get_cmap('BrownYellow'), vmin=1, vmax=100)
    ax.add_collection(lc)

    ax.axvline(27.5, color='w')
    ax.axvline(32.5, color='w')

    # Get moulin positions
    n_nodes = nodes.shape[0]
    moulin_xy_index = np.loadtxt('../../data/moulins/mesh_refinement/moulins_%05d.txt' % n_nodes).astype(int)
    moulin_xy = moulin_xy_index[:, 1:3]
    ax.plot(moulin_xy[:, 0]/1e3, moulin_xy[:, 1]/1e3, linestyle='', marker='x', color='k')


    ax.set_xlim([0, 100])
    ax.set_ylim([0, 25])
    ax.set_aspect('equal')



    out.close()
    fig.savefig(figname, dpi=600)

if __name__=='__main__':
    cases = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    for caseid in cases:
        fname = '../RUN/output_%03d_steady.nc' % caseid
        figname = 'steady_%03d.png' % caseid

        plot_steady(fname, figname)
    plt.show()
