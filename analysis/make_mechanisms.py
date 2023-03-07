"""

Call plot_mechanisms.py

"""

from matplotlib import pyplot as plt

from plot_mechanisms import plot_mechanisms
import defaults

## Case 00: Flat topo, SHMIP forcing
cases = [101, 102, 103, 104, 105]
fnames = ['../glads/00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = '00_mechanisms.png'
models = defaults.labels
fig_00 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.0137799761144865)


plt.show()
