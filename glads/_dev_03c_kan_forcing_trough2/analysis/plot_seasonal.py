import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib import gridspec

import cmocean

import GladsPlot as gplt
from palettes.code import palettes, tools

fname = '../RUN/output_001_seasonal.nc'
tslice = 365 + 190

cmap1 = palettes.get_cmap('BrownGray').reversed()
cmap2 = cmocean.cm.dense
cmap = tools.join_cmaps(cmap1, cmap2, average=0)

out = nc.Dataset(fname)

phi = out['phi'][:].data.T
N = out['N'][:].data.T

bed = np.vstack(out['bed'][:].data)
phi_bed = 9.81*1000*bed

pw = phi - phi_bed
ff = pw/(N + pw)

Q = np.abs(out['Q'][:].data.T)

nodes = out['nodes'][:].data.T/1e3
connect = out['connect'][:].data.T.astype(int) - 1
connect_edge = out['connect_edge'][:].data.T.astype(int) - 1

fig  = plt.figure(figsize=(8, 3))

mtri = Triangulation(nodes[:, 0], nodes[:, 1], connect)

gs = gridspec.GridSpec(2, 2, height_ratios=(5, 100),
    width_ratios=(100, 5), left=0.1, right=0.9,
    bottom=0.1, top=0.8, hspace=0)

ax = fig.add_subplot(gs[1, 0])
cax1 = fig.add_subplot(gs[0, 0])
cax2 = fig.add_subplot(gs[1, 1])

pc = ax.tripcolor(mtri, ff[:, tslice], cmap=cmap, vmin=-1, vmax=1)

lc = gplt.plot_edge_data(nodes, connect_edge, Q[:, tslice],
    palettes.get_cmap('BrownYellow'), vmin=1, vmax=100)
ax.add_collection(lc)
cb1 = fig.colorbar(lc, cax=cax1, orientation='horizontal')
cax1.xaxis.tick_top()
cax1.xaxis.set_label_position('top')

cb2 = fig.colorbar(pc, cax=cax2, orientation='vertical')

cb1.set_label(r'$Q~({\rm{m^3~s^{-1}}})$')
cb2.set_label(r'$p_{\rm{w}}/p_{\rm{i}}$')

Qticks = np.linspace(0, 100, 5)
Qticks[0] = 1
cb1.set_ticks(Qticks)

ax.set_aspect('equal')
ax.set_xlim([0, 100])
ax.set_ylim([0, 25])
ax.set_yticks([0, 12.5, 25])

plt.show()
