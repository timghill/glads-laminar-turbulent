"""
Nice 3D plot of bed and surface elevations
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import tri
from matplotlib.transforms import Bbox
from matplotlib.gridspec import GridSpec
import netCDF4 as nc
import cmocean
from palettes.code import palettes

dmesh = nc.Dataset('../glads/data/mesh/mesh_04.nc')
nodes = dmesh['tri/nodes'][:].data.T
connect = dmesh['tri/connect'][:].data.T.astype(int) - 1

out_trough = nc.Dataset('../glads/02_shmip_forcing_synth_topo/RUN/output_203_seasonal.nc')

out_valley = nc.Dataset('../glads/02a_shmip_forcing_valley_topo/RUN/output_003_seasonal.nc')


bed_flat = 350 + 0*nodes[:, 0]
bed_trough = out_trough['bed'][:].data
bed_valley = out_valley['bed'][:].data

surf = 6*(np.sqrt(nodes[:, 0] + 5e3) - np.sqrt(5e3)) + 390

mesh_tri = tri.Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)

bed_levels = np.arange(0, 600, 50) - 25
thick_levels = np.arange(0, 2100, 100)

## Plot bed elevation and ice thickness

fig2 = plt.figure(figsize=(7, 4.5))

gs = GridSpec(4, 2, height_ratios=(10, 100, 100, 100),
    left=0.1, bottom=0.15, right=0.97, top=0.8,
    hspace=0.05, wspace=0.075)

axs = np.array([[fig2.add_subplot(gs[i+1,j]) for j in range(2)] for i in range(3)]).T

bed_cmap = palettes.get_cmap('BrownEarth').reversed()
thick_cmap = palettes.get_cmap('BlueIce')

bed_pcolor = axs[0, 0].tricontourf(mesh_tri, bed_flat, cmap=bed_cmap, levels=bed_levels)
thick_pcolor = axs[1, 0].tricontourf(mesh_tri, surf - bed_flat, cmap=thick_cmap,  levels=thick_levels)

axs[0, 1].tricontourf(mesh_tri, bed_trough, cmap=bed_cmap, levels=bed_levels)
axs[1, 1].tricontourf(mesh_tri, surf - bed_trough, cmap=thick_cmap, levels=thick_levels)

axs[0, 2].tricontourf(mesh_tri, bed_valley, cmap=bed_cmap,  levels=bed_levels)
axs[1, 2].tricontourf(mesh_tri, surf - bed_valley, cmap=thick_cmap, levels=thick_levels)

for ax in axs.flat:
    ax.set_aspect('equal')
    ax.set_xlim([0, 100])
    ax.set_ylim([0, 25])
    ax.set_yticks([0, 12.5, 25])

axT = axs.T
axT[0, 0].set_xticklabels([])
axT[0, 1].set_xticklabels([])

axT[1, 0].set_xticklabels([])
axT[1, 1].set_xticklabels([])
# axT[0, 2].set_xticklabels([])

axT[0, 1].set_yticklabels([])
axT[1, 1].set_yticklabels([])
axT[2, 1].set_yticklabels([])

axT[0, 0].set_ylabel('y (km)')
axT[1, 0].set_ylabel('y (km)')
axT[2, 0].set_ylabel('y (km)')

axT[2, 0].set_xlabel('x (km)')
axT[2, 1].set_xlabel('x (km)')

cax1 = fig2.add_subplot(gs[0, 0])
cax2 = fig2.add_subplot(gs[0, 1])
cbar1 = fig2.colorbar(bed_pcolor, cax=cax1, orientation='horizontal')
cax1.xaxis.tick_top()
cax1.xaxis.set_label_position('top')

cbar2 = fig2.colorbar(thick_pcolor,cax=cax2, orientation='horizontal')
cax2.xaxis.tick_top()
cax2.xaxis.set_label_position('top')

# cticks = np.array(['%d'%z for z in np.arange(0, 2100, 100)])
# cticks[1::5] = ''
# cticks[2::5] = ''
# cticks[3::5] = ''
# cbar2.ax.set_xticklabels(cticks)

cbar1.set_label('Bed Elevation (m)')
cbar2.set_label('Thickness (m)')

fig2.savefig('bed_thickness.png', dpi=600)

plt.show()