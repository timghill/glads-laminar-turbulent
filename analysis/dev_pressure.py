"""

Call plot_pressure_maps_timeseries.py to make pressure figures

"""

import numpy as np
from matplotlib import pyplot as plt

from plot_pressure_maps_separate_timeseries import plot_pressure_maps_timeseries
import defaults

t_ticks = [1 + 4/12, 1 + 6/12, 1 + 8/12, 1 + 10/12]
# t_ticklabels = ['4', '6', '8', '10']
t_ticklabels = ['May', 'July', 'Sep', 'Nov']
t_lim = [t_ticks[0], t_ticks[-1]]
t_xlabel = 'Month'

KAN_opts = {'t_ticklabels':t_ticklabels,
            't_xlabel':t_xlabel,
            't_ticks':t_ticks,
            't_lim':t_lim}

synth_opts={'t_ticklabels':t_ticklabels[:-1],
            't_xlabel':t_xlabel,
            't_ticks':t_ticks[:-1],
            't_lim':[1 + 3/12, 1 + 9/12]}

# =============================================================================
## Synthetic forcing suite

## Case 00: Flat topo, synthetic forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/dev_00_synth_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_scaling.png'
fig_00 = plot_pressure_maps_timeseries(fnames, figname, melt_forcing='SHMIPadj', Qmin=1, Qmax=100,
    **synth_opts)

## Case 00: Flat topo, synthetic forcing
cases = [101, 102, 103, 104, 105]
fnames = ['../glads/dev_00_synth_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_10X.png'
fig_00 = plot_pressure_maps_timeseries(fnames, figname, melt_forcing='SHMIPadj', Qmin=1, Qmax=100,
    **synth_opts)


## Case 01: KAN forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/dev_01_kan_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_KAN.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, melt_forcing='KAN', Qmin=1, Qmax=100,
    **KAN_opts)

## Case 01: KAN forcing
cases = [101, 102, 103, 104, 105]
fnames = ['../glads/dev_01_kan_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_KAN_10X.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, melt_forcing='KAN', Qmin=1, Qmax=100,
    **KAN_opts)

## Case 01: KAN forcing
cases = [201, 202, 203, 204, 205]
fnames = ['../glads/dev_01_kan_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_KAN_20X.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, melt_forcing='KAN', Qmin=1, Qmax=100,
    **KAN_opts)

## Case 01: KAN forcing
cases = [301, 302, 303, 304, 305]
fnames = ['../glads/dev_01_kan_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_KAN_30X.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, melt_forcing='KAN', Qmin=1, Qmax=100,
    **KAN_opts)

## Case 01: KAN forcing
cases = [401, 402, 403, 404, 405]
fnames = ['../glads/dev_01_kan_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_KAN_40X.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, melt_forcing='KAN', Qmin=1, Qmax=100,
    **KAN_opts)

## Case 01: KAN forcing
cases = [501, 502, 503, 504, 505]
fnames = ['../glads/dev_01_kan_forcing_scaling/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/dev_scaling_KAN_50X.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, melt_forcing='KAN', Qmin=1, Qmax=100,
    **KAN_opts)


plt.show()
