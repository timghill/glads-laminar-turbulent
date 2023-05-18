"""

Simple plot of flux parameterizations

"""

import numpy as np
from matplotlib import pyplot as plt

from defaults import colors

# Physical constants
nu = 1.79e-6
rhow = 1000
rhoi = 910
g = 9.8

# Parameters
k_s = 0.1       # Laminar conductivity
omega = 1/2000  # Transition parameter
h_bed = 0.5

# Hydraulic potential scaling
phi_min = rhow*g*350 + rhoi*g*40
phi_max = rhow*g*350 + rhoi*g*1520
gradphi = (phi_max - phi_min)/100e3

# Conductivity scaling
h3 = nu/omega/k_s/gradphi
hcrit = h3**(1/3)
print('Transition h:', hcrit)
k_t = lambda alpha: k_s * h3**(1 - alpha/3) * gradphi**(2 - 3/2)

# Compute flux
h = np.logspace(-3, 0)
q_turb_54 = k_t(5/4)*h**(5/4)*gradphi**(0.5)
q_turb_32 = k_t(3/2)*h**(3/2)*gradphi**(0.5)

q_lam = k_s*h**3*gradphi

q_tran_54 = nu/omega * (h_bed/h)**(3 - 2*5/4) * ( -1 + np.sqrt(1 + 4*omega/nu * (h/h_bed)**(3-2*5/4)*k_s*h**3*gradphi))
q_tran_32 = nu/omega * (h_bed/h)**(3 - 2*3/2) * ( -1 + np.sqrt(1 + 4*omega/nu * (h/h_bed)**(3-2*3/2)*k_s*h**3*gradphi))

fig, ax = plt.subplots()

# ax.loglog(h/hcrit, omega/nu * q_turb_54, color=colors[0], linewidth=2, label='Turbulent 5/4')
# ax.loglog(h/hcrit, omega/nu * q_turb_32, color=colors[1], linewidth=1, label='Turbulent 3/2')
# ax.loglog(h/hcrit, omega/nu * q_lam, color=colors[2], linewidth=0.75, linestyle='dashed', zorder=5, label='Laminar')
# ax.loglog(h/hcrit, omega/nu * q_tran_54, color=colors[3], linewidth=2.5, label='Transition 5/4')
# ax.loglog(h/hcrit, omega/nu * q_tran_32, color=colors[4], linewidth=1.25, label='Transition 3/2')

ax.plot(h/hcrit, omega/nu * q_turb_54, color=colors[0], linewidth=2, label='Turbulent 5/4')
ax.plot(h/hcrit, omega/nu * q_turb_32, color=colors[1], linewidth=1, label='Turbulent 3/2')
ax.plot(h/hcrit, omega/nu * q_lam, color=colors[2], linewidth=0.75, linestyle='dashed', zorder=5, label='Laminar')
ax.plot(h/hcrit, omega/nu * q_tran_54, color=colors[3], linewidth=2.5, label='Transition 5/4')
ax.plot(h/hcrit, omega/nu * q_tran_32, color=colors[4], linewidth=1.25, label='Transition 3/2')
ax.legend()

fig.subplots_adjust(right=0.95)
fig.subplots_adjust(top=0.95)

# ax.set_ylim([1e-2, 1e2])
# ax.set_xlim([1e-1, 1e1])

ax.set_xlim([0, 2])
ax.set_ylim([0, 5])

ax.set_xlabel(r'$\tilde h$')
ax.set_ylabel(r'$\omega {\rm{Re}} = \frac{\omega q}{\nu}$')
ax.grid()

fig.savefig('flux_parameterizations_linear.png', dpi=600)

plt.show()

