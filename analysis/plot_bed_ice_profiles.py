"""

Plot bed and ice surface elevation profiles

"""

import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.transforms import Bbox

from mpl_toolkits.mplot3d.art3d import Poly3DCollection

#  class mpl_toolkits.mplot3d.art3d.Poly3DCollection(verts, *args, zsort='average', shade=False, lightsource=None, **kwargs)[source]

# Domain topography
bed_fun = lambda x: 350 + 0*x
surf_fun = lambda x: 6*(np.sqrt(x  +5e3) - np.sqrt(5e3)) + 390

# Graphics
# ice_color = 'lightcyan'
ice_color = 'lightblue'
bed_color = 'peru'
sky_color = 'whitesmoke'
line_color = 'k'
# band_color = 'lightslategrey'
band_color = 'steelblue'
figname = 'domain_profile.png'
figsize=(6, 6)

bands = [15, 30, 70]
band_width = 5
band_background = [0.5, 0.5, 0.5, 0.25]
# band_color = 'dimgray'

# Compute topo on x profile
x = np.linspace(0, 100e3, 101)
y = np.linspace(0, 25e3, 5)
[xx, yy] = np.meshgrid(x, y)
bed = bed_fun(xx)
surf = surf_fun(xx)

# Make the figure
fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(projection='3d', computed_zorder=False)

# Plot bed
ax.plot_surface(xx/1e3, yy/1e3, bed, color=bed_color, edgecolor='none', zorder=0)
xs = [0, 100]
ys = [0, 0]
zzs = 0*xs
[xxs, yys] = np.meshgrid(xs, ys)
print(xxs)
print(yys)
zzs = np.array([[0, 0], [350, 350]])
ax.plot_surface(xxs, yys, zzs, color=bed_color, zorder=1, antialiased=False, alpha=1)

xxs = np.array([[0, 0], [0, 0]])
yys = np.array([[0, 25], [0, 25]])
zzs = np.array([[0, 0], [350, 350]])
ax.plot_surface(xxs, yys, zzs, color=bed_color, zorder=1, antialiased=False, alpha=1)

# Plot surface
[xx1, yy1] = np.meshgrid(x, [0, 0])
z1 = surf_fun(xx1)
z1[0] = 390
print(xx1)
print(yy1)
print(z1)
ax.plot_surface(xx/1e3, yy/1e3, surf, color=ice_color, edgecolor=ice_color, alpha=1, zorder=0, antialiased=False, linewidth=0.25)
ax.plot_surface(xx1/1e3, yy1/1e3, z1, color=ice_color, edgecolor=ice_color, alpha=1, zorder=0, antialiased=False, linewidth=0.25)

# Plot bands
for xb in bands:
    xs = np.linspace(xb-2.5, xb+2.5, 101)
    ys = np.array([0, 25])
    [xxs, yys] = np.meshgrid(xs, ys)

    zs = np.array([[350, surf_fun(xb*1e3)], [350, surf_fun(xb*1e3)]])

    ax.plot([xb, xb], [0, 25], [surf_fun(xb*1e3), surf_fun(xb*1e3)], color=line_color, zorder=5, linewidth=2)
    ax.plot([xb, xb], [0, 0], [390, surf_fun(xb*1e3)], color=line_color, zorder=5, linewidth=2)
    zzs = surf_fun(xxs*1e3)
    ax.plot_surface(xxs, yys, zzs+1, zorder=2, color=band_color, alpha=1, antialiased=False)
    z3 = zzs.copy()
    z3[0] = 350
    ax.plot_surface(xxs, yys*0, z3, zorder=2, color=band_color, alpha=1, antialiased=False)

# Plot moulins
# dmesh = nc.Dataset('../glads/data/mesh/mesh_04.nc')
moulins = np.loadtxt('../glads/data/moulins/moulins_068.txt')
print(moulins)
moulinx = moulins[:, 1]
mouliny = moulins[:, 2]
moulinz = surf_fun(moulinx)

ax.plot(moulinx/1e3, mouliny/1e3, moulinz, marker='.', color='k', zorder=6, linestyle='')


ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_yticks([0, 12.5, 25])
ax.set_zticks([0, 500, 1000, 1500, 2000])

ax.set_xlabel('Distance from terminus (km)', labelpad=20)
ax.set_ylabel('Distance (km)')
ax.zaxis.set_rotate_label(False)  # disable automatic rotation
ax.set_zlabel('z (m)', rotation=90)

ax.view_init(elev=24, azim=-125) #Works!
ax.set_box_aspect((4, 1, 1))
ax.set_aspect('equalxy')


# ax.set_box_aspect(1)
# ax.set_aspect('equal')

# ax.set_aspect([10, 10, 10])


# fig, ax = plt.subplots(figsize=figsize, projection='3d')

# ax.fill_between(x/1e3, bed, 0, facecolor=bed_color)
# ax.fill_between(x/1e3, surf, bed, facecolor=ice_color, edgecolor=line_color)
# ax.axhline(350, color='dodgerblue', linewidth=2)

# for xb in bands:
#     xx_band = np.linspace(xb - band_width/2, xb + band_width/2, 21)
#     yy_band = surf_fun(xx_band*1e3)
#     ax.fill_between(xx_band, yy_band, 350, facecolor=band_background)
#     ax.plot([xb, xb], [350, surf_fun(xb*1e3)], color=band_color)


# ax.set_xlim([0, 100]); ax.set_ylim([0, xb00])
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# ax.set_xlabel('Distance from terminus')
# ax.set_ylabel('Elevation (m)')

# fig.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.23)
# ax.set_position(Bbox.from_bounds(-0.6, 0.5, 2, 0.65))
ax.set_position(Bbox.from_extents(0.1, 0.4, 1, 1.15))

ax2 = fig.add_subplot()
ax2.set_position(Bbox.from_bounds(0.15, 0.1, 0.8, 0.35))

ax2.set_xlabel('Elevation (m)')
ax2.set_ylabel(r'Density (km$^{-2}$)')
ax2.grid()

fig.savefig(figname, dpi=600)

plt.show()

