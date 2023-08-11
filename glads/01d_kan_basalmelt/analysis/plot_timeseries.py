import numpy as np
import netCDF4 as nc

from matplotlib import pyplot as plt

fname = '../RUN/output_001_seasonal.nc'

xbands = [15, 30, 70]

with nc.Dataset(fname, 'r') as out:
    phi = out['phi'][:].data.T
    N = out['N'][:].data.T

    bed = np.vstack(out['bed'][:].data)
    phi_bed = 9.81*1000*bed

    pw = phi - phi_bed
    ff = pw/(N + pw)
    nodes = out['nodes'][:].data.T/1e3
    time = out['time'][:].data

tt = (time/86400/365 - 101)*12

fig, ax = plt.subplots()
for i,xb in enumerate(xbands):
    node_mask = np.logical_and(nodes[:, 0]<=(xb + 2.5), nodes[:, 0]>=(xb - 2.5))
    ff_band = np.mean(ff[node_mask, :], axis=0)
    ax.plot(tt, ff_band, label='x = %d km' % xb)

ax.grid()
ax.legend()
ax.set_xlabel('Month')
ax.set_ylabel(r'$p_{\rm{w}}/p_{\rm{i}}$')
ax.set_xlim([3, 9])
plt.show()
