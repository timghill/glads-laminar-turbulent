"""
postproc.py

Plot floatation fraction, Q
Plot distributed h_s, Q

All for just one timestep

"""

import os

import numpy as np
import netCDF4 as nc

from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.gridspec import GridSpec
import cmocean

import GladsPlot as gplt
from palettes.code import palettes, tools

def do_postproc(filename, test=False):

    out = nc.Dataset(filename)
    path, fname = os.path.split(filename)
    figname = fname.replace('.nc', '.png')
    
    time = out['time'][:].data
    tt = out['time'][:].data - 100*365*86400
    N = out['N'][:].data.T
    phi = out['phi'][:].data.T

    ix = np.arange(len(tt))
    taxis = tt/86400

    rho_i = out['para/rho_i'][:].data; rho_w = out['para/rho_w'][:].data
    g = out['para/g_grav'][:].data

    thick = np.vstack(out['thick'][:].T)
    p_ice = rho_w*g*thick
    p_w = p_ice - N

    f = p_w/p_ice

    hs = out['h_sheet'][:].data.T

    Q = np.abs(out['Q'][:].data.T)

    h_bed = float(out['para/h_bed'][:].data)


    dt_days = int((tt[1] - tt[0])/86400)
    tii = 180//dt_days
    fig = plt.figure(figsize=(8, 5))
    gs = GridSpec(3, 2, height_ratios=(10, 100, 100), width_ratios=[100, 4],
        hspace=0.15, left=0.05, bottom=0.05, wspace=0.05, top=0.9)

    axs = np.array([fig.add_subplot(gs[i+1, 0]) for i in range(2)])
    caxs = np.array([fig.add_subplot(gs[j+1, 1]) for j in range(2)])

    nodes = out['nodes'][:].data.T
    connect = out['connect'][:].data.T - 1
    connect_edge = out['connect_edge'][:].data.T.astype(int) - 1
    mtri = Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)
    titles = ['S', 'Q']
    lc_vmin = 10
    lc_vmax = 100
    
    ax1 = axs[0]
    ax2 = axs[1]

    fcolor = ax1.tripcolor(mtri, f[:, tii], vmin=0, vmax=2, cmap=cmocean.cm.curl)
    ax1.set_aspect('equal')
    ax1.set_xlim([0, 100])
    ax1.set_ylim([0, 25])

    lc = gplt.plot_edge_data(nodes/1e3, connect_edge, Q[:, tii], palettes.get_cmap('BrownYellow'),
        vmin=lc_vmin, vmax=lc_vmax)
    ax1.add_collection(lc)

    
    hcolor = ax2.tripcolor(mtri, hs[:, tii], vmin=0, vmax=h_bed, cmap=cmocean.cm.haline)
    ax2.tricontour(mtri, hs[:, tii], levels=[out['para/h_bed'][:].data], colors='k')
    ax2.set_aspect('equal')
    ax2.set_xlim([0, 100])
    ax2.set_ylim([0, 25])

    lc = gplt.plot_edge_data(nodes/1e3, connect_edge, Q[:, tii], palettes.get_cmap('BrownYellow'),
        vmin=lc_vmin, vmax=lc_vmax)
    ax2.add_collection(lc)
    
    cb_f = fig.colorbar(fcolor, cax=caxs[0])
    cb_h = fig.colorbar(hcolor, cax=caxs[1])
    cax_top = fig.add_subplot(gs[0, 0])
    cb_top = fig.colorbar(lc, cax=cax_top, orientation='horizontal')

    cb_f.set_label(r'$p_{\rm{w}}/p_{\rm{i}}$')
    cb_h.set_label(r'$h_{\rm{s}}~(\rm{m})$')
    cb_top.set_label(r'$Q~(\rm{m}^3~\rm{s}^{-1})$')
    
    cax_top.xaxis.set_label_position('top')
    cax_top.xaxis.tick_top()

    ax1.set_xticklabels([])

    # ax1.text(0.5, 1.5, 'Seasonal sensitivity %s' % fname, fontweight='bold', transform=ax1.transAxes,
    #     fontsize=12, ha='center')

    if not test:
        fig.savefig(figname, dpi=600)



if __name__ == '__main__':
    fname_base = '../RUN/output_%03d_seasonal.nc'
    for jobid in [501, 502, 503]:
        fname = fname_base % jobid
        do_postproc(fname)
