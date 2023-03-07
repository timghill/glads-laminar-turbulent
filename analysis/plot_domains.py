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


bed_flat = 350 + 0*nodes[:, 0]
bed_trough = out_trough['bed'][:].data

surf = 6*(np.sqrt(nodes[:, 0] + 5e3) - np.sqrt(5e3)) + 390

mesh_tri = tri.Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)

fig = plt.figure(figsize=(7, 4))
ax1 = fig.add_subplot(1, 2, 1, projection='3d', box_aspect=(100, 25, 1.6*20), elev=20, azim=-120)
ax2 = fig.add_subplot(1, 2, 2, projection='3d', box_aspect=(100, 25, 1.6*20), elev=20, azim=-120)
psurf = ax1.plot_trisurf(mesh_tri, surf, cmap=palettes.get_cmap('BlueIce'), vmin=0, vmax=2000, edgecolors=(0, 0, 0, 0.5), linewidths=0.2)
ax1.plot_trisurf(mesh_tri, bed_flat, cmap=palettes.get_cmap('BrownEarth'), vmin=0, vmax=400)
ax1.set_xlim([0, 100])
ax1.set_ylim([0, 25])
ax1.set_zlim([0, 1950])
ax1.set_yticks([0, 12.5, 25])

ax1.set_xlabel('x (km)')
ax1.set_ylabel('y (km)')

cbar = fig.colorbar(psurf, location='bottom', aspect=40, ax=ax1)
cbar.set_label('Surface elevation (m)')

# ax1.set_position(Bbox.from_bounds(-0.95, -0.1, 1.7, 1.2))
# ax1.set_position(Bbox.from_bounds(-0.45, 0.05, 2, 1.2))


ax2.plot_trisurf(mesh_tri, surf, cmap=palettes.get_cmap('BlueIce'), vmin=0, vmax=2000, edgecolors=(0, 0, 0, 0.5), linewidths=0.2)
bsurf = ax2.plot_trisurf(mesh_tri, bed_trough, cmap=palettes.get_cmap('BrownEarth'), vmin=0, vmax=400)
ax2.set_xlim([0, 100])
ax2.set_ylim([0, 25])
ax2.set_zlim([0, 1950])

ax2.set_xlabel('x (km)')
ax2.set_ylabel('y (km)')
ax2.set_yticks([0, 12.5, 25])
# ax2.set_position(Bbox.from_bounds(-0.45, 0.05, 2, 1.2))

cbar = fig.colorbar(bsurf, location='bottom', aspect=40, ax=ax2)
cbar.set_label('Bed elevation (m)')

# plt.tight_layout()
fig.savefig('shmip_domain')


## Plot bed elevation and ice thickness

fig2 = plt.figure(figsize=(7, 3))

gs = GridSpec(3, 2, height_ratios=(10, 100, 100),
    left=0.1, bottom=0.15, right=0.97, top=0.8,
    hspace=0.05, wspace=0.075)

axs = np.array([[fig2.add_subplot(gs[i+1,j]) for j in range(2)] for i in range(2)])

bed_cmap = palettes.get_cmap('BrownEarth').reversed()
thick_cmap = palettes.get_cmap('BlueIce')

bed_pcolor = axs[0, 0].tricontourf(mesh_tri, bed_flat, cmap=bed_cmap, vmin=0, vmax=400, levels=np.arange(0, 450, 50))

thick_pcolor = axs[1, 0].tricontourf(mesh_tri, surf - bed_flat, cmap=thick_cmap, vmin=0, vmax=2000, levels=np.arange(0, 2100, 100))

axs[0, 1].tricontourf(mesh_tri, bed_trough, cmap=bed_cmap, vmin=0, vmax=400, levels=np.arange(0, 450, 50))

axs[1, 1].tricontourf(mesh_tri, surf - bed_trough, cmap=thick_cmap, vmin=0, vmax=2000, levels=np.arange(0, 2100, 100))

for ax in axs.flat:
    ax.set_aspect('equal')
    ax.set_xlim([0, 100])
    ax.set_ylim([0, 25])
    ax.set_yticks([0, 12.5, 25])

axs[0, 0].set_xticklabels([])
axs[0, 1].set_xticklabels([])

axs[0, 1].set_yticklabels([])
axs[1, 1].set_yticklabels([])

axs[0, 0].set_ylabel('y (km)')
axs[1, 0].set_ylabel('y (km)')

axs[1, 0].set_xlabel('x (km)')
axs[1, 1].set_xlabel('x (km)')

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