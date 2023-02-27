import numpy as np
from matplotlib import pyplot as plt
import netCDF4 as nc

# files = ['../RUN/output_501_steady.nc', '../RUN/output_502_steady.nc', '../RUN/output_503_steady.nc']
files = ['../RUN/output_501_seasonal.nc', '../RUN/output_502_seasonal.nc', '../RUN/output_503_seasonal.nc']

fig, ax = plt.subplots()
labels = ['Turbulent', 'Laminar', 'Transition']
colors = ['darkmagenta', 'cornflowerblue', 'mediumslateblue']
for i, f in enumerate(files):
    D = nc.Dataset(f)
    phi = D['phi'][:].data.T
    pw = phi
    N = D['N'][:].data.T
    f = pw/(N + pw)
    nodes = D['nodes'][0].data/1e3    
    ax.scatter(nodes, f[:, 180+365], 10, alpha=0.5, label=labels[i], color=colors[i])
    ax.scatter(nodes, f[:, -1], 10, alpha=0.5, color=colors[i])

ax.legend(markerscale=3)
ax.set_ylim([0, 1])
ax.set_xlim([0, 100])
ax.grid()

fig.savefig('profiles_50x_seasonal.png', dpi=600)
plt.show()

