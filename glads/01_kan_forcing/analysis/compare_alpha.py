"""

Compare alpha >>3 with alpha=3, turbulent, transition models

"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib import gridspec
import scipy.interpolate
import netCDF4 as nc
import cmocean
import GladsPlot as gplt
from palettes.code import palettes

cases = [2, 3, 6]
colors = np.array([[0.579, 0.677, 0.781, 1],
                   [0.250, 0.250, 0.250, 1],
                   [0.500, 0.500, 0.500, 1]])
labels = ['Turbulent 3/2', 'Laminar 3', 'Laminar 4']
band = 30e3
band_width = 5e3
tslice = 365 + 190
# tslice = 365 + 168
Qmin = 1
Qmax = 100

ts_fig, ts_ax = plt.subplots()
map_fig = plt.figure(figsize=(7, 6))
map_gs = gridspec.GridSpec(4, 2, height_ratios=(10, 100, 100, 100),
    width_ratios=(100, 3), hspace=0.1, wspace=0.1, left=0.1)
map_axs = np.array([map_fig.add_subplot(map_gs[i+1, 0]) for i in range(len(cases))])
map_caxs = np.array([map_fig.add_subplot(map_gs[0, 0]),
                     map_fig.add_subplot(map_gs[1:, 1])])

fig2, ax2 = plt.subplots()

for i,caseid in enumerate(cases):
    fname = '../RUN/output_%03d_seasonal.nc' % caseid
    print(fname)
    with nc.Dataset(fname, 'r') as out:
        phi = out['phi'][:].data.T
        N = out['N'][:].data.T

        bed = np.vstack(out['bed'][:].data)
        Q = np.abs(out['Q'][:].data.T)
        
        connect = out['connect'][:].data.T.astype(int) - 1
        nodes = out['nodes'][:].data.T
        elements = out['elements'][:].data.T
        connect_edge = out['connect_edge'][:].data.T.astype(int) - 1
        times = out['time'][:].data.T/86400/365 - 101
        qxy = out['qs'][:].data.T
        qs = np.sqrt(qxy[:, 0]**2 + qxy[:, 1]**2)
        hs = out['h_sheet'][:].data.T
        ks = out['para/cond_s'][:].data.T
    
    phi_bed = 1000*9.8*bed
    pw = phi - phi_bed
    ff = pw/(N + pw)

    node_mask = np.logical_and(nodes[:, 0]>(band - band_width/2),
                                nodes[:, 0]<=(band + band_width/2))
    ff_mean = np.mean(ff[node_mask, :], axis=0)
    ts_ax.plot(times*12, ff_mean, color=colors[i], label=labels[i])

    mtri = Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)
    map_ax = map_axs[i]
    pc = map_ax.tripcolor(mtri, ff[:, tslice], cmap=cmocean.cm.dense, vmin=0, vmax=1)
    map_ax.set_aspect('equal')
    map_ax.set_xlim([0, 100])
    map_ax.set_ylim([0, 25])
    map_ax.set_yticks([0, 12.5, 25])
    map_ax.set_xticks([])

    map_ax.tricontour(mtri, ff[:, tslice], colors='w', levels=[1], linewidths=0.5)

    map_ax.text(0.95, 0.1, labels[i], ha='right', color='w', fontweight='bold',
        transform=map_ax.transAxes)

    lc = gplt.plot_edge_data(nodes/1e3, connect_edge, Q[:, tslice],
            palettes.get_cmap('BrownYellow'), vmin=Qmin, vmax=Qmax)
    map_ax.add_collection(lc)

    interpo = scipy.interpolate.LinearNDInterpolator(nodes, hs)
    h_el = interpo(elements)
    if caseid==2:
        gradphi = (qs/ks/h_el**(3./2.))**2
    elif caseid==3:
        gradphi = qs/ks/h_el**3
    elif caseid==6:
        gradphi = qs/ks/h_el**4
    
    k_eff = (qs/h_el**(5./4.)/gradphi**0.5)

    el_mask = np.logical_and(elements[:, 0]<=(band + band_width/2),
                             elements[:, 0]>(band - band_width/2))
    k_eff_mean = np.nanmean(k_eff[el_mask, :], axis=0)
    ax2.plot(times*12, k_eff_mean, color=colors[i], label=labels[i])

ts_ax.grid(linestyle=':', linewidth=0.5)
ts_ax.set_xlabel('Month')
ts_ax.set_ylabel('Flotation fraction')
ts_ax.set_xlim([4, 11])
ts_ax.legend()
ts_ax.axvline((tslice - 365)/365 * 12, color='k', linewidth=0.5)

ax2.grid(linestyle=':', linewidth=0.5)
ax2.set_xlabel('Month')
ax2.set_ylabel(r'$k_{\rm{eff}}$')
ax2.set_xlim([4, 11])
ax2.legend()

cbar_top = map_fig.colorbar(lc, cax=map_caxs[0], orientation='horizontal')
cbar_right = map_fig.colorbar(pc, cax=map_caxs[1])

map_caxs[0].xaxis.tick_top()
map_caxs[0].xaxis.set_label_position('top')
map_axs[-1].set_xlabel('Distance from terminus (km)')
map_axs[1].set_ylabel('Distance across-glacier (km)')

map_axs[-1].set_xticks([0, 20, 40, 60, 80, 100])

cbar_top.set_label(r'$Q~({\rm{m}}^3~{\rm{s}}^{-1})$')
cbar_right.set_label('Flotation fraction')

ts_fig.subplots_adjust(top=0.95, right=0.95)
fig2.subplots_adjust(top=0.95, right=0.95)

map_fig.savefig('laminar_alpha_4_map.png', dpi=600)
ts_fig.savefig('laminar_alpha_4_ts.png', dpi=600)
fig2.savefig('laminar_alpha_4_cond.png', dpi=600)

plt.show()
