"""
Plot distributed melt defined by KAN_L 2014 timeseries
"""

import datetime as dt

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

import cmocean

# Parameter
lapse = -0.005
DDF = 0.01

AWS = np.loadtxt('KAN_L_2014_temp_clipped.txt', delimiter=',')
days = AWS[:, 0]
sl_temp = AWS[:, 1]

tt = np.array([dt.datetime(2014, 1, 1) + dt.timedelta(day) for day in days])

xx = np.linspace(0, 100e3, 101)
z = 6*(np.sqrt(xx + 5000) - np.sqrt(5000)) + 390

temp_dist = np.zeros((len(z), len(tt)))
for i in range(len(tt)):
    temp_dist[:, i] = sl_temp[i] + lapse*z

melt_dist = temp_dist*DDF
melt_dist[melt_dist<0] = 0

fig = plt.figure()
gs = GridSpec(2, 2, width_ratios=(100, 3), hspace=0.05, wspace=0.05,
    right=0.875, top=0.9, left=0.1)

[XX, YY] = np.meshgrid(xx, days)

ax = fig.add_subplot(gs[0, 0])
pc = ax.pcolor(XX/1e3, YY, melt_dist.T, cmap=cmocean.cm.amp)
cax = fig.add_subplot(gs[0, 1])
cbar = fig.colorbar(pc, cax=cax, orientation='vertical')
cbar.set_label('Melt (m w.e. day$^{-1}$)')

# ax.contour(melt_dist.T, levels=[0], colors='k')
ax.set_xticklabels([])
ax.set_xlim([0, 100])

ax.set_ylabel('Day of year')
ax.set_title(r'Lapse Rate = %.4f$^{\circ}\rm{C}~\rm{m}^{-1}$' % lapse)

ax2 = fig.add_subplot(gs[1, 0])
ax2.plot(xx/1e3, np.sum(melt_dist.T, axis=0))
ax2.grid()
ax2.set_xlabel('x (km)')
ax2.set_xlim([0, 100])
ax2.set_ylabel('Cumulative melt (m w.e.)')

fig.savefig('KAN_dist_melt.png', dpi=600)

plt.show()
