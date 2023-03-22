"""

Plot Moody diagram using Colebrook equation

"""

import numpy as np

from scipy.optimize import newton
from scipy import interpolate

from matplotlib import pyplot as plt
import netCDF4 as nc

pattern = '../glads/00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'
cases = [201, 202, 203, 204, 205]
models = ['Turbulent 5/4', 'Turbulent 3/2', 'Laminar', 'Transition 5/4', 'Transition 3/2']
fnames = [pattern%id for id in cases]
print(fnames)

colors = np.array([[0.420, 0.510, 0.620, 1],
                   [0.579, 0.677, 0.781, 1],
                   [0.500, 0.500, 0.500, 1],
                   [0.859, 0.683, 0.275, 1],
                   [0.929, 0.835, 0.408, 1]])

# PARAMETERS
omega = 1/2000
nu = 1.796e-6
k_lam = 0.1
h_0 = 0.5
rhow = 1000

fig, ax = plt.subplots()

for ii in range(len(fnames)):
    out = nc.Dataset(fnames[ii], 'r')

    # Read parameters
    k = out['para/cond_s'][:].data

    h = out['h_sheet'][:].data.T
    qxy = out['qs'][:].data.T
    qs = np.sqrt(qxy[:, 0]**2 + qxy[:, 1]**2)
    Re = qs/nu
    nodes = out['nodes'][:].data.T
    elements = out['elements'][:].data.T
    h_node = out['h_sheet'][:].data.T
    interpolator = interpolate.LinearNDInterpolator(nodes, h_node)
    hs = interpolator(elements)
    
    if models[ii]=='Turbulent 5/4':
        gradphi = (qs/k/hs**(5/4))**2

        h_bed = float(out['para/h_bed'][:].data)
        max_f_D = h_bed**0.5 / rhow / k**2
    elif models[ii]=='Turbulent 3/2':
        gradphi = (qs/k/hs**(3/2))**2
    elif models[ii]=='Laminar':
        gradphi = (qs/k/hs**(3))**1
    elif models[ii]=='Transition 5/4':
        h_bed = float(out['para/h_bed'][:].data)
        gradphi = ((qs + omega/nu * (hs/h_bed)**(1/2) * qs**2)/k/hs**(3))**1
    elif models[ii]=='Transition 3/2':
        gradphi = ((qs + omega/nu * qs**2)/k/hs**(3))**1

    f_D = hs**3*gradphi/rhow/Re**2/nu**2

    step = 5
    f_D = f_D[::step, ::step].flatten()
    Re = Re[::step, ::step].flatten()

    ax.scatter(Re.flatten(), f_D.flatten(), color=colors[ii], label=models[ii], s=1, alpha=0.3)

ax.set_xscale('log')
ax.set_yscale('log')

ax.grid(linestyle=':')

ax.set_xlabel(r'${\rm{Re}} = \frac{q}{\nu}$', fontsize=12)
ax.set_ylabel(r'$f_{\rm{D}} = \frac{h^3 |\nabla \phi|}{\rho_{\rm{w}} \nu^2 {\rm{Re}}^2}$', fontsize=12)
ax.legend(markerscale=5, loc='upper right', framealpha=1, edgecolor='none')
ax.axvline(1/omega, color='k', linewidth=1)
# ax.axhline(max_f_D, color='k', linewidth=1)
plt.tight_layout()

fig.savefig('moody_scatter.png', dpi=600)

# plt.show()
