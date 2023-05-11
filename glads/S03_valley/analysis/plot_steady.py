import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
import netCDF4 as nc

out = nc.Dataset('../RUN/output_001_steady.nc')

connect = out['connect'][:].data.T - 1
nodes = out['nodes'][:].data.T
mtri = Triangulation(nodes[:, 0], nodes[:, 1], connect)

fig, ax = plt.subplots()
pc = ax.tripcolor(mtri, out['h_sheet'][:].data.T[:, -1])
fig.colorbar(pc)
ax.set_aspect('equal')
ax.set_title('h_sheet')


fig, ax = plt.subplots()
pc = ax.tripcolor(mtri, out['N'][:].data.T[:, -1])
fig.colorbar(pc)
ax.set_aspect('equal')
ax.set_title('N')


plt.show()
