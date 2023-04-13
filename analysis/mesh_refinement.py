"""

Plot results of mesh refinement tests

"""

import numpy as np
import netCDF4 as nc

from matplotlib import pyplot as plt

def plot_mesh_refinement(fnames, figname, figsize=(6, 4)):
    """Some string
    """

    fig, ax1 = plt.subplots(figsize=figsize)

    n_meshes = len(fnames)
    n_nodes = np.zeros(n_meshes)
    mean_edge_lengths = np.zeros(n_meshes)
    ff_arr = np.zeros(n_meshes)
    Q_arr = np.zeros(n_meshes)
    clock_arr = np.zeros(n_meshes)
    for (i,fname) in enumerate(fnames):
        # Read model output
        print(fname)
        out = nc.Dataset(fname, 'r')


        # Read mesh file
        dmesh_fname = '../glads/data/mesh/mesh_refinement_%02d.nc' % (i+1)
        dmesh = nc.Dataset(dmesh_fname, 'r')
        nodes = dmesh['tri/nodes'][:].data.T
        connect = dmesh['tri/connect'][:].data.T.astype(int) - 1
        connect_edge = dmesh['tri/connect_edge'][:].data.T.astype(int) - 1
        area_nodes = dmesh['tri/area_nodes'][:].data.T

        # Compute averaged quantities
        xb = 30e3       # Set x position for averaging
        bw = 5e3        # Band width for averaging
        # node_mask = np.abs(nodes[:, 0] - xb)<=bw/2
        node_mask = nodes[:, 0]>=0

        # Floatation fraction: check pressure steady state
        phi = out['phi'][-1].data
        rhow = 1000
        g = 9.8
        phi0 = rhow*g*out['bed'][:].data
        pw = phi - phi0

        N = out['N'][-1].data
        pi = N + pw
        ff = pw/pi
        ff = ff[node_mask]
        ff_mean = np.sum(ff*area_nodes)/np.sum(area_nodes)

        # Check flux or discharge steady state
        Q = out['Q'][-1].data
        nodex1 = nodes[connect_edge[:, 0], 0]
        nodex2 = nodes[connect_edge[:, 1], 0]
        channel_mask_pos = np.logical_and(nodex1>=xb, nodex2<xb)
        channel_mask_neg = np.logical_and(nodex2>=xb, nodex1<xb)

        Q_pos = Q[channel_mask_pos]
        Q_neg = Q[channel_mask_neg]
        Qtot = np.sum(Q_pos) - np.sum(Q_neg)

        # Walltime
        clock = out['model/wallclock'][:].data

        n_nodes[i] = nodes.shape[0]
        mean_edge_lengths[i] = n_nodes[i]**2
        ff_arr[i] = ff_mean
        Q_arr[i] = Qtot
        clock_arr[i] = clock
    
    h1, = ax1.semilogx(n_nodes, ff_arr, marker='o', label=r'$p_{\rm{w}}/p_{\rm{i}}$', color='r')
    ax1.grid()
    ax1.set_xlabel('Nodes')
    
    ax2 = ax1.twinx()
    h2, = ax2.semilogx(n_nodes, Q_arr, marker='^', label='Q (m$^3$ s$^{-1}$)', color='b')
    ax2.yaxis.tick_left()
    ax2.yaxis.set_label_position('left')
    ax2.spines.left.set_position(('axes', -0.1))
    
    ax3 = ax1.twinx()
    h3, = ax3.semilogx(n_nodes, clock_arr/3600, marker='x', label='Time (h)', color='k')

    ax1.legend(handles=(h1, h2, h3))

    plt.show()


        


if __name__=='__main__':
    cases = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    fnames = ['../glads/S00_mesh_refinement/RUN/output_%03d_steady.nc' % caseid for caseid in cases]
    figname = 'S00_mesh_refinement.png'
    plot_mesh_refinement(fnames, figname)