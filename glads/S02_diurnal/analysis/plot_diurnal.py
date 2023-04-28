"""
Plot diurnal simulations
"""

import numpy as np
from matplotlib import pyplot as plt
import netCDF4 as nc

def _get_ff(ncout):
    phi = ncout['phi'][:].data.T
    rhow = 1000
    g = 9.8
    bed = np.vstack(ncout['bed'][:].data)
    phi_elev = rhow*g*bed
    pw = phi - phi_elev

    N = ncout['N'][:].data.T
    ff = pw/(N + pw)

    return ff

def diurnal(cases, figname):
    fig, ax = plt.subplots(figsize=(8, 4))

    colors = np.array([[0.579, 0.677, 0.781, 1],
                    [0.199, 0.328, 0.492, 1],
                    [0.500, 0.500, 0.500, 1],
                    [0.929, 0.835, 0.408, 1],
                    [0.859, 0.683, 0.275, 1]])
    labels = ['Turbulent 5/4', 'Turbulent 3/2', 'Laminar', 'Transition 5/4', 'Transition 3/2']

    for j in range(len(cases)):
        fname_steady = '../RUN/output_%03d_steady.nc' % cases[j]
        fname_diurnal = '../RUN/output_%03d_seasonal.nc' % cases[j]
        out_steady = nc.Dataset(fname_steady, 'r')
        out_diurnal = nc.Dataset(fname_diurnal, 'r')

        ff_steady = _get_ff(out_steady)
        ff_diurnal = _get_ff(out_diurnal)

        nodes = out_diurnal['nodes'][:].data.T
        node_mask = np.abs(nodes[:, 0] - 30e3)<=2.5e3

        tt_steady = out_steady['time'][:].data
        tt_diurnal = out_diurnal['time'][:].data

        tt_steady = (tt_steady - tt_steady[-1])/86400
        tt_diurnal = (tt_diurnal - tt_diurnal[0])/86400
        ax.plot(tt_steady, np.mean(ff_steady[node_mask], axis=0) - np.mean(ff_steady[node_mask], axis=0)[-1], color=colors[j])
        ax.plot(tt_diurnal, np.mean(ff_diurnal[node_mask], axis=0) -np.mean(ff_steady[node_mask], axis=0)[-1], color=colors[j], label=labels[j])

    ax.grid(linestyle=':', linewidth=0.5)
    ax.set_xlim([-0.5, 5.5])
    ax.set_ylim([-0.15, 0.15])
    ax.legend(bbox_to_anchor=(-0.1, 1.02, 1.2, 0.1), mode='expand', ncol=5, frameon=False)
    ax.set_xlabel('Days')
    ax.set_ylabel(r'$\Delta p_{\rm{w}}/p_{\rm{i}}$')

    fig.savefig(figname, dpi=600)

if __name__=='__main__':
    cases = [101, 102, 103, 104, 105]

    figname = 'diurnals.png'

    diurnal(cases, figname)
    plt.show()
