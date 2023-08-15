"""

Plot KAN sea level temperature

"""

import numpy as np
from matplotlib import pyplot as plt


# Graphics
figname = 'KAN_forcing.png'
figsize=(6, 3)
color = 'k'

# Get melt data
KAN_day_temp = np.loadtxt('../../glads/data/kan_l_melt/KAN_L_2014_temp_clipped.txt', delimiter=',')

# Time axis for plotting
tt_years = np.arange(3/12, 12/12, 1/365/2)
tt_days = tt_years*365

KAN_tt = KAN_day_temp[:, 0]
KAN_temp = KAN_day_temp[:, 1]

temp_interp = np.interp(tt_days, KAN_tt, KAN_temp, left=np.nan, right=np.nan)

fig, ax = plt.subplots(figsize=figsize)
temp_interp[tt_years>(10/12)] = np.nan
ax.plot(tt_years*12 + 1, temp_interp, color=color)
# ax.set_xlim([5, 10])
ax.grid(linestyle=':', linewidth=0.5)

month_labels = ['May', 'June', 'July', 'Aug', 'Sep', 'Oct', 'Nov']
ax.set_xticks([5, 6, 7, 8, 9, 10, 11])
ax.set_xlim([4.75, 11.25])
ax.set_xticklabels(month_labels)
ax.set_ylabel('Temperature (C)')

ax.axhline(0, color='k', linewidth=0.5)

fig.subplots_adjust(bottom=0.1, top=0.95, right=0.98, left=0.125)

fig.savefig(figname, dpi=600)

plt.show()