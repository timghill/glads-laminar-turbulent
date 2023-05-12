import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
import netCDF4 as nc

import cmocean

import GladsPlot as gplt

case = 605
out = nc.Dataset('../RUN/output_%d_steady.nc' % (case))

connect = out['connect'][:].data.T - 1
nodes = out['nodes'][:].data.T

bed = np.vstack(out['bed'][:].data.T)
thick = np.vstack(out['thick'][:].data.T)
rhow = 1000
g = 9.8
nu = 1.79e-6
phi_z = rhow*g*bed

phi = out['phi'][:].data.T
N = out['N'][:].data.T
pw = phi - phi_z
pi = 910*g*thick
ff = pw/pi

Q = np.abs(out['Q'][:].data.T)

qxy = out['qs'][:].data.T
qs = np.sqrt(qxy[:, 0]**2 + qxy[:, 1]**2)
elements = out['elements'][:].data.T

h_sheet = out['h_sheet'][:].data.T

time = out['time'][:].data.T
tt = time/86400/365

Pa_min = 0
Pa_max = 3e6

Qmin = 1e-1
Qmax = 10


mtri = Triangulation(nodes[:, 0], nodes[:, 1], connect)

fig, ax = plt.subplots()
pc = ax.tripcolor(mtri, out['h_sheet'][:].data.T[:, -1])
fig.colorbar(pc)
ax.set_aspect('equal')
ax.set_title('h_sheet')


fig, ax = plt.subplots(figsize=(6, 2.5))
pc = ax.tripcolor(mtri, (phi - phi_z)[:, -1], vmin=-Pa_max, vmax=Pa_max, cmap=cmocean.cm.balance)
fig.colorbar(pc)
ax.set_aspect('equal')
ax.set_title('p_w')

line_cmap = cmocean.cm.matter
connect_edge = out['connect_edge'][:].data.T.astype(int) - 1
lc = gplt.plot_edge_data(nodes, connect_edge, Q[:, -1],
            line_cmap, vmin=Qmin, vmax=Qmax)
ax.add_collection(lc)

fig.savefig('SHMIP_valley_pw_%d.png' % case, dpi=600)

fig, ax = plt.subplots()
pc = ax.tripcolor(mtri, N[:, -1], vmin=Pa_min, vmax=Pa_max)
fig.colorbar(pc)
ax.set_aspect('equal')
ax.set_title('N')

fig, ax = plt.subplots()
pc = ax.tripcolor(mtri, qs[:, -1]/nu, vmin=0, vmax=4000, cmap=cmocean.cm.curl)
fig.colorbar(pc)
ax.set_aspect('equal')
ax.set_title('Re')

fig, ax = plt.subplots(figsize=(6, 2.5))
pc = ax.tripcolor(mtri, pi.flatten(), vmin=-Pa_max, vmax=Pa_max, cmap=cmocean.cm.balance)
fig.colorbar(pc)
ax.set_aspect('equal')
ax.set_title('p_i')
fig.savefig('SHMIP_valley_overburden_%d.png' % case, dpi=600)

fig, ax = plt.subplots()
ax.plot(tt, np.mean(h_sheet, axis=0))
ax.set_title('h_sheet')

fig, ax = plt.subplots()
# ax.plot(tt, np.mean(ff, axis=0))
# ax.set_title('Floatation fraction')
ax.plot(tt, np.mean(phi - phi_z, axis=0))
ax.plot(tt, np.mean(pi - N, axis=0))


plt.show()
