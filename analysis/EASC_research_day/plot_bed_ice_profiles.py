"""

Plot bed and ice surface elevation profiles

"""

import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D

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
line_color = 'steelblue'
figname = 'domain_profile.png'
figsize=(6, 3)

bands = [20]
band_width = 5
band_background = [0.5, 0.5, 0.5, 0.25]
band_color = 'dimgray'

# Compute topo on x profile
x = np.linspace(0, 100e3, 101)
y = np.linspace(0, 25e3, 5)
[xx, yy] = np.meshgrid(x, y)
bed = bed_fun(xx)
surf = surf_fun(xx)

# Make the figure
fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(projection='3d')

ax.plot_surface(xx/1e3, yy/1e3, bed, color=bed_color, edgecolor='none')

# xs = np.array([20, 20, 20, 20])
# ys = np.array([0, 25, 25, 0])
xs = np.linspace(20-2.5, 20+2.5, 101)
ys = np.array([0, 25])
[xxs, yys] = np.meshgrid(xs, ys)

zs = np.array([[350, surf_fun(20e3)], [350, surf_fun(20e3)]])

# zs = np.array([350, 350, surf_fun(20e3), surf_fun(20e3)])

# ax.plot_surface(xxs, yys, zs, color='k')

# verts = np.array([xx.flatten(), yy.flatten(), surf.flatten()])
# verts = np.zeros((len(xx.flatten()), 3))
# verts[:, 0] = xx.flatten()/1e3
# verts[:, 1] = yy.flatten()/1e3
# verts[:, 2] = surf.flatten()
# p3col = Poly3DCollection([verts])
# ax.add_collection(p3col)
ax.plot_surface(xx/1e3, yy/1e3, surf, color=ice_color, edgecolor='none', alpha=0.8)

[xx1, yy1] = np.meshgrid(x, [0, 0])
z1 = surf_fun(xx1)
z1[0] = 390
print(xx1)
print(yy1)
print(z1)
ax.plot_surface(xx1/1e3, yy1/1e3, z1, color=ice_color, edgecolor='none', alpha=0.8)

ax.plot([20, 20], [0, 25], [surf_fun(20e3) + 50, surf_fun(20e3) + 50], color=line_color, zorder=5)
ax.plot([20, 20], [0, 0], [390, surf_fun(20e3)], color=line_color, zorder=5)
zzs = surf_fun(xxs*1e3)
ax.plot_surface(xxs, yys, zzs+20, zorder=10, color=line_color, alpha=0.5)
z3 = zzs.copy()
z3[0] = 350
ax.plot_surface(xxs, yys*0 - 0.01, z3, zorder=10)
# verts = zip(xx)



ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_yticks([0, 12.5, 25])
ax.set_zticks([0, 500, 1000, 1500, 2000])

ax.set_xlabel('Distance from terminus (km)')
ax.set_ylabel('Distance (km)')
ax.set_zlabel('Elevation (m)')

ax.view_init(elev=40, azim=-135) #Works!


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


# ax.set_xlim([0, 100]); ax.set_ylim([0, 2000])
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# ax.set_xlabel('Distance from terminus')
# ax.set_ylabel('Elevation (m)')

# fig.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.23)

fig.savefig(figname, dpi=600)

plt.show()

