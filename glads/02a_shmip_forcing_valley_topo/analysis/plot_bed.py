import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
import netCDF4 as nc

out = nc.Dataset('../RUN/output_303_seasonal.nc')
nodes = out['nodes'][:].data.T
connect = out['connect'][:].data.T - 1
fig, ax = plt.subplots()

mtri = Triangulation(nodes[:, 0], nodes[:, 1], connect)

bed = ax.tripcolor(mtri, out['bed'][:].data)
fig.colorbar(bed)

ax.set_aspect('equal')

plt.show()
