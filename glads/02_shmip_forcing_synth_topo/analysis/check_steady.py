"""
Check if seasonal simulations have reached a quasi-steady-state
"""

import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt

fname_pattern = '../RUN/output_%03d_seasonal.nc'
cases = [921, 924, 922, 923, 925]
n_cases = len(cases)
figname = 'convergence.png'

rhow = 1000
g = 9.81

tindices = [365, 365*2]

dmesh = nc.Dataset('../../data/mesh/mesh_04.nc')
node_area = np.vstack(dmesh['tri/area_nodes'][:].data)
edge_length = np.vstack(dmesh['tri/edge_length'][:].data)


for i in range(n_cases):
    out = nc.Dataset(fname_pattern % cases[i])
    
    phi = out['phi'][:].data.T
    phi_m = rhow*g*np.vstack(out['bed'][:].data)
    p_w = phi - phi_m
    N = out['N'][:].data.T
    ff = p_w/(N + p_w);

    S = out['S_channel'][:].data.T

    time = out['time'][:]/86400/365 - 100

    out.close()

    print('Case %s:' % (fname_pattern % cases[i]))
    print('Mean, Abs Max difference flotation fraction')
    mean_ff_1 = np.sum(node_area*ff[:, tindices[0]])/np.sum(node_area)
    mean_ff_2 = np.sum(node_area*ff[:, tindices[1]])/np.sum(node_area)
    delta_ff = mean_ff_2 - mean_ff_2
    
    print(delta_ff)
    
    max_delta = np.max(np.abs(ff[:, tindices[1]] - ff[:, tindices[0]]))
    print(max_delta)

    print('Mean, Abs Max difference channel area')
    delta_S = S[:, tindices[1]] - S[:, tindices[0]]
    print(np.mean(delta_S))
    print(np.max(np.abs(delta_S)))
    
    fig, ax = plt.subplots()
    mean_ff = np.sum(node_area*ff/np.sum(node_area), axis=0)
    time_mod = time % 1
    dt_mod = np.zeros(len(time_mod))
    dt_mod[1:] = time_mod[1:] - time_mod[:-1]
    mean_ff[dt_mod<0] = np.nan

    mean_S = np.sum(edge_length*S/np.sum(edge_length), axis=0)
    mean_S[dt_mod<0] = np.nan

    ax.plot(time % 1, mean_ff, label=r'$p_{\rm{w}}/p_{\rm{i}}$')
    ax.plot(time % 1, mean_S, label='$S$')
    ax.grid()
    ax.set_xlabel('Years from transition')
    ax.set_title(fname_pattern % cases[i])

    ax.legend()


plt.show()

    


