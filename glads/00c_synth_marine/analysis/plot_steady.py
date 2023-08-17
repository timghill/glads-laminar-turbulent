import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib import gridspec

import cmocean

import GladsPlot as gplt
from palettes.code import palettes

out = nc.Dataset('../RUN/output_001_steady.nc')

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

fig  = plt.figure()

mtri = Triangulation(nodes[:, 0], nodes[:, 1], connect)

gs = gridspec.GridSpec(2, 2, height_ratios=(5, 100),
    width_ratios=(100, 5))

ax = fig.add_subplot(gs[1, 0])
cax1 = fig.add_subplot(gs[0, 0])
cax2 = fig.add_subplot(gs[1, 1])

pc = ax.tripcolor(mtri, ff[:, -1], cmap=cmocean.cm.dense, vmin=0, vmax=1)

lc = gplt.plot_edge_data(nodes, connect_edge, Q[:, -1],
    palettes.get_cmap('BrownYellow'), vmin=1, vmax=100)
ax.add_collection(lc)
cb1 = fig.colorbar(lc, cax=cax1, orientation='horizontal')
cb2 = fig.colorbar(pc, cax=cax2, orientation='vertical')

ax.set_aspect('equal')
ax.set_xlim([0, 100])
ax.set_ylim([0, 25])

plt.show()
