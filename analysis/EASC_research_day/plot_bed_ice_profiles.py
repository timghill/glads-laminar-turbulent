"""

Plot bed and ice surface elevation profiles

"""

import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors

# Domain topography
bed_fun = lambda x: 350 + 0*x
surf_fun = lambda x: 6*(np.sqrt(x  +5e3) - np.sqrt(5e3)) + 390

# Graphics
ice_color = 'lightcyan'
bed_color = 'peru'
sky_color = 'whitesmoke'
line_color = 'lightblue'
figname = 'domain_profile.png'
figsize=(6, 2)

bands = [20]
band_width = 5
band_background = [0.5, 0.5, 0.5, 0.25]
band_color = 'dimgray'

# Compute topo on x profile
x = np.linspace(0, 100e3, 101)
bed = bed_fun(x)
surf = surf_fun(x)

# Make the figure
fig, ax = plt.subplots(figsize=figsize)

ax.fill_between(x/1e3, bed, 0, facecolor=bed_color)
ax.fill_between(x/1e3, surf, bed, facecolor=ice_color, edgecolor=line_color)
ax.axhline(350, color='dodgerblue', linewidth=2)

for xb in bands:
    xx_band = np.linspace(xb - band_width/2, xb + band_width/2, 21)
    yy_band = surf_fun(xx_band*1e3)
    ax.fill_between(xx_band, yy_band, 350, facecolor=band_background)
    ax.plot([xb, xb], [350, surf_fun(xb*1e3)], color=band_color)


ax.set_xlim([0, 100]); ax.set_ylim([0, 2000])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_xlabel('Distance from terminus')
ax.set_ylabel('Elevation (m)')

fig.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.23)

fig.savefig(figname, dpi=600)

plt.show()

