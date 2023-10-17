import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import netCDF4 as nc

import defaults

def plot_issm(mat_fnames, issm_fnames, figname,
    colors=defaults.colors
    ):

    # x_bands = [15e3, 30e3, 70e3]
    xb = 30e3

    # Set up the figure and axes
    fig = plt.figure(figsize=(7, 4))
    gs = GridSpec(2, 3, left=0.1, right=0.98, bottom=0.12, top=0.95)
    axs = np.array([[fig.add_subplot(gs[i,j]) for j in range(3)] for i in range(2)])
    alphabet = ['a', 'b', 'c', 'd', 'e']
    for ii in range(len(mat_fnames)):
        mfname = mat_fnames[ii]
        ifname = issm_fnames[ii]

        with nc.Dataset(mfname, 'r') as mdata:
            nodes = mdata['nodes'][:].data.T
            phi = mdata['phi'][:].data.T
            N = mdata['N'][:].data.T
            phi_0 = 9.8*1000*np.vstack(mdata['bed'][:].data)
            pw = phi - phi_0
            ff = pw/(N + pw)
            mtt = mdata['time'][:].data/86400/365 - 101
        
        with nc.Dataset(ifname, 'r') as idata:
            N_issm = idata['N'][:].data.T
            phi_issm = idata['phi'][:].data.T
            pw_issm = phi_issm - np.vstack(idata['phi_bed'][:].data)
            ff_issm = pw_issm/(N_issm + pw_issm)
            itt = idata['time'][:].data.T

        itt = itt - 9

        ax = axs.flat[ii]
        nodemask = np.logical_and(
            nodes[:,0]<=(xb+2.5e3), nodes[:,0]>=(xb-2.5e3))
        # nodemask = nodes[:,0]>=0
        ax.plot(mtt, np.mean(ff[nodemask, :], axis=0),
            color=colors[ii])
        ax.plot(itt, np.mean(ff_issm[nodemask, :], axis=0),
            color=colors[ii], linestyle=':')
        
        ax.set_xlim([0, 1])
        ax.set_ylim([0.4, 1])
        ax.text(0.05, 0.95, alphabet[ii], va='top', fontweight='bold')
        ax.grid(linewidth=0.5, linestyle=':')
    
    axs.flat[-1].set_visible(False)

    for ax in axs[0]:
        ax.set_xticklabels([])
    
    for ax in axs[:, 1:].flat:
        ax.set_yticklabels([])
    

    fig.text(0.01, 0.5, r'$p_{\rm{w}}/p_{\rm{i}}$', va='center', rotation='vertical')
    fig.text(0.5, 0.02, 'Years', va='bottom')
    axs.flat[4].legend(labels=('Matlab', 'ISSM'),
        bbox_to_anchor=(1.05, 0, 1, 1), mode='expand',
        frameon=False)

    if figname:
        fig.savefig(figname, dpi=400)
        

if __name__ == '__main__':
    mat_fnames = ['/home/tghill/scratch/laminar-turbulent/glads/00_synth_forcing/RUN/output_%03d_seasonal.nc'%i for i in range(1,6)]
    issm_fnames = ['/home/tghill/scratch/laminar-turbulent/issm/00_synth_forcing/RUN/output_%03d.nc'%i for i in range(1,6)]
    print(mat_fnames)
    print(issm_fnames)
    figname = 'figures/aux/issm_comparison.png'
    plot_issm(mat_fnames, issm_fnames, figname)

    plt.show()
