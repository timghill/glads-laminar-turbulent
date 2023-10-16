"""
Plot Reynolds number profiles
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import netCDF4 as nc

import cmocean
from palettes.code import palettes, tools

figsize=(7, 3.5)

gs_kwargs=dict(wspace=0.075, hspace=0.2, 
        left=0.09, right=0.95, bottom=0.12,
        top=0.95)

labels = ['Turbulent 5/4', 'Turbulent 3/2', 'Laminar',
            'Transition 5/4', 'Transition 3/2']

def main(fnames, figname, figsize=figsize,
    gs_kwargs=gs_kwargs, map_cmap=None, Re_ylim=[0, 4e3],
    labels=labels
    ):
    n_cases = len(fnames)
    # Set up figure
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(3, 2, **gs_kwargs)
    cax_gs = gridspec.GridSpecFromSubplotSpec(4, 1,
        subplot_spec=gs[1, 1], height_ratios=(100, 30, 30, 100),
        hspace=0.2)
    

    dmesh = nc.Dataset('/home/tghill/projects/def-gflowers/tghill/laminar-turbulent/glads/data/mesh/mesh_04.nc')
    nodes = dmesh['tri/nodes'][:].data.T
    connect = dmesh['tri/connect'][:].data.T
    # elements = dmesh['tri/elements'][:].data.T
    area = dmesh['tri/area'][:].data.T
    edges = dmesh['tri/edge_midpoints'][:].data.T


    if map_cmap is None:
        cm1 = palettes.get_cmap('BrownGray')
        cm2 = palettes.get_cmap('blue-8').reversed()
        z1 = int(200*2000/Re_ylim[1])
        z2 = 200 - z1
        map_cmap = tools.join_cmaps(cm1, cm2, N1=z1, N2=z2, average=10)
    
    ii_rows = [0, 0, 1, 2, 2]
    ii_cols = [0, 1, 0, 0, 1]
    axs = np.array([fig.add_subplot(gs[i,j]) for i,j in zip(ii_rows, ii_cols)])
    alphabet = ['a', 'b', 'c', 'd', 'e']

    # Set up width-averaging
    xe = np.linspace(0, 100e3, 101)
    xc = 0.5*(xe[1:] + xe[:-1])

    alphabet = ['a', 'b', 'c', 'd', 'e']

    # Start reading the data
    for ii in range(n_cases):
        fname = fnames[ii]
        print(fname)
        ax = axs[ii]

        out = nc.Dataset(fname)
        elements = out['elements'][:].data.T

        Q = np.abs(out['S_channel'][:].data.T)

        q = out['qs'][:].data.T
        qs = np.sqrt(q[:, 0, :]**2 + q[:, 1, :]**2)
        Re = qs/float(out['para/nu'][:].data)
        time = out['time'][:].data
        time = (time - time[0])/365/86400
        x = xc/1e3

        [tt, xx] = np.meshgrid(time, x)

        # We need to width-average Re...
        Re_avg = np.zeros((xc.shape[0], len(time)))
        for jj in range(len(xc)):
            node_mask = np.logical_and(elements[:, 0]>=xe[jj], elements[:, 0]<=xe[jj+1])
            Re_avg[jj, :] = np.mean(Re[node_mask, :], axis=0)
        
        xrep = np.zeros(Q.shape)
        xrep[:, :] = np.vstack(edges[:, 0])
        xrep[Q<1] = np.nan
        Qlimit = np.nanmax(xrep, axis=0)
        print(Re_avg.shape)
        pc = ax.pcolormesh(tt, xx, Re_avg, cmap=map_cmap, vmin=Re_ylim[0], vmax=Re_ylim[1], shading='nearest')
        ax.set_xlim([1 + 4/12, 1 + 10/12])

        ax.text(0.05, 0.95, alphabet[ii], fontweight='bold', transform=ax.transAxes,
            va='top', ha='left')
        ax.text(0.95, 0.95, labels[ii], transform=ax.transAxes,
            va='top', ha='right')

    cax1 = fig.add_subplot(cax_gs[1])
    cb = fig.colorbar(pc, cax=cax1, orientation='horizontal')
    cb.set_label(r'Re')

    for ax in axs:
        ax.set_xticks([1+4/12, 1+6/12, 1+8/12, 1+10/12])
        ax.set_xticklabels(['May', 'July', 'Sep', 'Nov'])

    for ax in axs[:3]:
        ax.set_xticklabels([])
    
    for ax in axs[[1,4]]:
        ax.set_yticklabels([])

    axs[2].set_ylabel('Distance from terminus (km)')
    axs[-1].set_xticks([1+6/12, 1+8/12, 1+10/12])
    axs[-1].set_xticklabels(['July', 'Sep', 'Nov'])
    
    for ax in axs[-2:]:
        ax.set_xlabel('Month')

    fig.savefig(figname, dpi=300)

        

if __name__=='__main__':

    cases = [1, 2, 3, 4, 5]
    fpattern = '/home/tghill/scratch/laminar-turbulent/glads/00_synth_forcing/RUN/output_%03d_seasonal.nc'
    fnames = [fpattern%caseid for caseid in cases]
    figname = 'figures/main/00_synth_Re.png'
    main(fnames, figname, Re_ylim=(0, 4e3))

    cases = [1, 2, 3, 4, 5]
    fpattern = '/home/tghill/scratch/laminar-turbulent/glads/00a_shmip_forcing/RUN/output_%03d_seasonal.nc'
    fnames = [fpattern%caseid for caseid in cases]
    figname = 'figures/aux/00a_shmip_Re.png'
    main(fnames, figname, Re_ylim=(0, 4e3))

    cases = [1, 2, 3, 4, 5]
    fpattern = '/home/tghill/scratch/laminar-turbulent/glads/01_kan_forcing/RUN/output_%03d_seasonal.nc'
    fnames = [fpattern%caseid for caseid in cases]
    figname = 'figures/supplement/01_KAN_Re.png'
    main(fnames, figname, Re_ylim=(0, 4e3))

    cases = [1, 2, 3, 4, 5]
    fpattern = '/home/tghill/scratch/laminar-turbulent/glads/01a_kan_adj_forcing/RUN/output_%03d_seasonal.nc'
    fnames = [fpattern%caseid for caseid in cases]
    figname = 'figures/aux/01a_KAN_adj_Re.png'
    main(fnames, figname, Re_ylim=(0, 4e3))

    cases = [1, 2, 3, 4, 5]
    fpattern = '/home/tghill/scratch/laminar-turbulent/glads/S01a_parameter_sensitivity/RUN/output_%03d_seasonal.nc'
    fnames = [fpattern%caseid for caseid in cases]
    figname = 'figures/aux/S01a_Re.png'
    main(fnames, figname, Re_ylim=(0, 4e3))

    cases = [1, 2, 3, 4, 5]
    fpattern = '/home/tghill/scratch/laminar-turbulent/glads/S01b_parameter_sensitivity/RUN/output_%03d_seasonal.nc'
    fnames = [fpattern%caseid for caseid in cases]
    figname = 'figures/aux/S01b_Re.png'
    main(fnames, figname, Re_ylim=(0, 4e3))

    
    plt.show()
