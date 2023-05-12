import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
import netCDF4 as nc

out = nc.Dataset('../RUN/output_005_steady.nc')

connect = out['connect'][:].data.T - 1
nodes = out['nodes'][:].data.T

bed = np.vstack(out['bed'][:].data.T)
thick = np.vstack(out['thick'][:].data.T)
print(bed.shape)
print(thick.shape)
rhow = 1000
g = 9.8
phi_z = rhow*g*bed

phi = out['phi'][:].data.T
N = out['N'][:].data.T
pw = phi - phi_z
pi = 910*g*thick
ff = pw/pi

h_sheet = out['h_sheet'][:].data.T

time = out['time'][:].data.T
tt = time/86400/365


mtri = Triangulation(nodes[:, 0], nodes[:, 1], connect)

fig, ax = plt.subplots()
pc = ax.tripcolor(mtri, out['h_sheet'][:].data.T[:, -1])
fig.colorbar(pc)
ax.set_aspect('equal')
ax.set_title('h_sheet')


fig, ax = plt.subplots()
pc = ax.tripcolor(mtri, (phi - phi_z)[:, -1])
fig.colorbar(pc)
ax.set_aspect('equal')
ax.set_title('phi stuff')

fig, ax = plt.subplots()
ax.plot(tt, np.mean(h_sheet, axis=0))
ax.set_title('h_sheet')

fig, ax = plt.subplots()
# ax.plot(tt, np.mean(ff, axis=0))
# ax.set_title('Floatation fraction')
ax.plot(tt, np.mean(phi - phi_z, axis=0))
ax.plot(tt, np.mean(pi - N, axis=0))


plt.show()
