"""

Call plot_Re to plot Reynolds number maps

"""

import numpy as np
from matplotlib import pyplot as plt

from plot_Re import plot_Re

## Case 00: Flat topo, synthetic forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/dev_00_synth_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_Re_day225.png'
fig_00 = plot_Re(fnames, figname, Qmin=0.01, Qmax=10, Re_ylim=(0, 4e3), tslice=365+225)

## Case 00: Flat topo, synthetic forcing
cases = [101, 102, 103, 104, 105]
fnames = ['../glads/dev_00_synth_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_Re_10X_day225.png'
fig_00 = plot_Re(fnames, figname, Qmin=0.01, Qmax=10, Re_ylim=(0, 4e3), tslice=365+225)



## Case 01: Flat topo, KAN_L forcing
cases = [1, 2, 3, 4, 5]
pattern = '../glads/dev_01_kan_forcing_scaling/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_KAN_Re.png'
fig_01 = plot_Re(fnames, figname, Qmin=1, Qmax=100, Re_ylim=(0, 4e3), tslice=365+174)

## Case 01: Flat topo, KAN_L forcing
cases = [401, 402, 403, 404, 405]
pattern = '../glads/dev_01_kan_forcing_scaling/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_KAN_Re_40X.png'
fig_01 = plot_Re(fnames, figname, Qmin=1, Qmax=100, Re_ylim=(0, 4e3), tslice=365+174)

## Case 01: Flat topo, KAN_L forcing
cases = [501, 502, 503, 504, 505]
pattern = '../glads/dev_01_kan_forcing_scaling/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_KAN_Re_50X.png'
fig_01 = plot_Re(fnames, figname, Qmin=1, Qmax=100, Re_ylim=(0, 4e3), tslice=365+174)


plt.show()
