"""
Compute water budget for GlaDS simulations
"""

import numpy as np
import netCDF4 as nc

from matplotlib import pyplot as plt

glads_fname = '../../glads/00_synth_forcing/RUN/output_001_seasonal.nc'
mesh_fname = '../../glads/data/mesh/mesh_04.nc'

# Read outputs
with nc.Dataset(glads_fname, 'r') as out:
    phi = out['phi'][:].data
    rhow = out['para/rho_w'][:].data
    g = out['para/g_grav'][:].data
    bed = out['bed'][:].data
    phi_bed = rhow*g*bed
    pw = phi - phi_bed

    hs = out['h_sheet'][:].data
    S = out['S_channel'][:].data
    e_v = out['para/e_v'][:].data
    time = out['time'][:].data

# Read mesh
with nc.Dataset(mesh_fname, 'r') as dmesh:
    node_area = dmesh['tri/area_nodes'][:].data
    element_area = dmesh['tri/area'][:].data
    edge_length = dmesh['tri/edge_length'][:].data

# Compute water budget terms
sheet_volume = np.sum(hs*node_area, axis=1)
channel_volume = np.sum(S*edge_length, axis=1)
storage_volume = np.sum(e_v*pw*node_area/rhow/g, axis=1)

total_volume = sheet_volume + channel_volume + storage_volume


fig, ax = plt.subplots()
tt = (time-101*365*86400)/365/86400
ax.semilogy(tt, sheet_volume)
ax.semilogy(tt, channel_volume)
ax.semilogy(tt, storage_volume)
# ax.plot(tt, total_volume)
ax.grid()
ax.legend(['Sheet', 'Channel', 'Storage'])
ax.set_xlim([0, 1])
ax.set_xlabel('Year')
ax.set_ylabel(r'Water volume (m$^3$)')

# Compute rates
d_dt = lambda x: (x[1:] - x[:-1])/(time[1:] - time[:-1])

d_mass = lambda x: (x - x[365])

fig, ax = plt.subplots()
# ax.plot(tt, d_mass(sheet_volume))
# ax.plot(tt, d_mass(channel_volume))
# ax.plot(tt, d_mass(storage_volume))
# ax.plot(tt, d_mass(total_volume))

ax.plot(tt[1:], d_dt(sheet_volume))
ax.plot(tt[1:], d_dt(channel_volume))
ax.plot(tt[1:], d_dt(storage_volume))
ax.plot(tt[1:], d_dt(total_volume))
ax.set_xlim([0, 1])
ax.grid()



plt.show()