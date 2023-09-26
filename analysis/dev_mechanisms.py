"""

Call plot_mechanisms.py

"""

from matplotlib import pyplot as plt

from plot_mechanisms import plot_mechanisms
import defaults

## Case 00a: Flat topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/dev_00_synth_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_00_mechanisms.png'
models = defaults.labels
fig_00 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.00954298166213413, tslice=530)

## Case 00a: Flat topo, SHMIP forcing
cases = [101, 102, 103, 104, 105]
fnames = ['../glads/dev_00_synth_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_00_mechanisms_10X.png'
models = defaults.labels
fig_00 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.00714915847943998, tslice=530)


## Case 1: KAN forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/dev_01_kan_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_01_mechanisms.png'
models = defaults.labels
fig_01 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.00954298166213413, tslice=365+174, labels=models)

## Case 1: KAN forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/dev_01_kan_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_01_mechanisms_40x.png'
models = defaults.labels
fig_01 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.0071, tslice=365+174, labels=models)

plt.show()
