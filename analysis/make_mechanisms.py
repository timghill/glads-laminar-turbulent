"""

Call plot_mechanisms.py

"""

from matplotlib import pyplot as plt

from plot_mechanisms import plot_mechanisms
import defaults

## Case 00a: Flat topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
prefix = '/home/tghill/projects/def-gflowers/tghill/laminar-turbulent/'
fnames = [prefix + 'glads/00_synth_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/main/00_mechanisms.png'
models = defaults.labels
fig_00 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.00714915847943998, tslice=530)

plt.show()

