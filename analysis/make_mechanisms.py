"""

Call plot_mechanisms.py

"""

from matplotlib import pyplot as plt

from plot_mechanisms import plot_mechanisms
import defaults

## Case 00a: Flat topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/00_synth_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = '00_mechanisms.png'
models = defaults.labels
fig_00 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.00954298166213413, tslice=530)

"""
## Case 00a: Flat topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/00a_shimp_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = '00a_mechanisms_shmip_forcing.png'
models = defaults.labels
fig_00 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.00954298166213413, tslice=530)
"""

plt.show()

