"""

Call plot_pressure_maps_timeseries.py to make pressure figures

"""

import numpy as np
from matplotlib import pyplot as plt

from plot_pressure_maps_separate_timeseries import plot_pressure_maps_timeseries
import defaults

t_ticks = [1 + 4/12, 1 + 6/12, 1 + 8/12, 1 + 10/12]
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

prefix = '/home/tghill/projects/def-gflowers/tghill/laminar-turbulent/'

# =============================================================================
## Synthetic forcing suite

## Case 00: Flat topo, synthetic forcing
cases = [1, 2, 3, 4, 5]
fnames = [prefix + 'glads/00_synth_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/main/00_pressure_seasonal.png'
fig_00 = plot_pressure_maps_timeseries(fnames, figname, melt_forcing='SHMIPadj', Qmin=1, Qmax=100,
    **synth_opts)

## Case 00a: Flat topo, standard SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = [prefix + 'glads/00a_shmip_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/supplement/00a_pressure_seasonal_shmip_forcing.png'
fig_00a = plot_pressure_maps_timeseries(fnames, figname, Qmin=10, Qmax=200, melt_forcing='SHMIP',
    **synth_opts)

## Case 00c: synthetic forcing, marine outlet
cases = [1, 2, 3, 4, 5]
fnames = [prefix + 'glads/00c_synth_marine/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/00c_pressure_seasonal_marine.png'
fig_00 = plot_pressure_maps_timeseries(fnames, figname, melt_forcing='SHMIPadj', Qmin=1, Qmax=100,
    **synth_opts)

# =============================================================================
## KAN forcing suite

## Case 01: Flat topo, KAN_L forcing
cases = [1, 2, 3, 4, 5]
pattern = prefix + 'glads/01_kan_forcing/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/main/01_pressure_seasonal.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, Qmax=100, melt_forcing='KAN',
     **KAN_opts, ff_ylim=[0, 1.75], ff_yticks=[0, 0.5, 1, 1.5])

## Case 01a: KAN increased forcing
cases = [1, 2, 3, 4, 5]
pattern = prefix + 'glads/01a_kan_adj_forcing/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname ='figures/supplement/01a_pressure_seasonal_KANadj.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, Qmin=10, Qmax=200, melt_forcing='KANadj',
     **KAN_opts,
     ff_ylim=[0, 2.75], ff_yticks=[0, 0.5, 1, 1.5, 2, 2.5])

## Case 01b: Flat topo, KAN_L forcing, reduced e_v
cases = [1, 2, 3, 4, 5]
pattern = prefix + 'glads/01b_kan_forcing_ev/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/supplement/01b_pressure_seasonal_ev.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, Qmax=100, melt_forcing='KAN',
     **KAN_opts,
     ff_ylim=[0, 2.75], ff_yticks=[0, 0.5, 1, 1.5, 2, 2.5])

## Case 01d: Basal melt rate
cases = [1, 2, 3, 4, 5]
pattern = prefix + 'glads/01d_kan_basalmelt/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/01d_pressure_seasonal_basalmelt.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, Qmax=100, melt_forcing='KAN',
     **KAN_opts, ff_ylim=[0, 1.75], ff_yticks=[0, 0.5, 1, 1.5])

# =============================================================================
## Bed topo suite

## Case 03a: Trough topo, KAN forcing
cases = [1, 2, 3, 4, 5]
pattern = prefix + 'glads/03a_kan_forcing_trough/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/supplement/03a_pressure_seasonal_trough.png'
fig_03 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, melt_forcing='KAN',
    **KAN_opts, ff_ylim=[0,1.75], ff_yticks=[0, 0.5, 1, 1.5])

## Case 03b: Valley topo, KAN forcing
cases = [1, 2, 3, 4, 5]
pattern = prefix + 'glads/03b_kan_forcing_valley/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/supplement/03b_pressure_seasonal_valley.png'
fig_03 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, melt_forcing='KAN',
    **KAN_opts, ff_ylim=[0,1.75], ff_yticks=[0, 0.5, 1, 1.5])

## Case 03c: Sinusoidal trough topo, KAN forcing
cases = [1, 2, 3, 4, 5]
pattern = prefix + 'glads/03c_kan_forcing_trough2/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/supplement/03c_pressure_seasonal_trough2.png'
fig_03c = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, Qmax=200, melt_forcing='KAN',
    **KAN_opts, ff_ylim=[0,1.75], ff_yticks=[0, 0.5, 1, 1.5])

# =============================================================================
## Parameter sensitivity
## Case S01: Parameter sensitivity

cases = [1, 2, 3, 4, 5]
pattern = prefix + 'glads/S01a_parameter_sensitivity/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/supplement/S01a_pressure_seasonal_params.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, Qmax=100, melt_forcing='KAN',
     **KAN_opts, ff_ylim=[0, 1.75], ff_yticks=[0, 0.5, 1, 1.5])


## Case S01: Parameter sensitivity
cases = [1, 2, 3, 4, 5]
pattern = prefix + 'glads/S01b_parameter_sensitivity/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/supplement/S01b_pressure_seasonal_params.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, Qmax=100, melt_forcing='KAN',
     **KAN_opts, ff_ylim=[0, 1.75], ff_yticks=[0, 0.5, 1, 1.5])


plt.show()
