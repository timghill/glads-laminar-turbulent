"""

Call plot_pressure_maps_timeseries.py to make pressure figures

"""

from matplotlib import pyplot as plt

from plot_pressure_maps_separate_timeseries import plot_pressure_maps_timeseries

t_ticks = [1 + 4/12, 1 + 6/12, 1 + 8/12, 1 + 10/12]
t_ticklabels = ['4', '6', '8', '10']
t_lim = [t_ticks[0], t_ticks[-1]]
t_xlabel = 'Month'

"""
## Case 00: Flat topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = '00_pressure_seasonal.png'
fig_00 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, melt_forcing='SHMIP',
    t_ticklabels=t_ticklabels[:-1], t_xlabel=t_xlabel, t_ticks=t_ticks[:-1], t_lim=[1 + 3/12, 1 + 9/12])
"""

## Case 00.2: Flat topo, SHMIP forcing
cases = [201, 202, 203, 204, 205]
fnames = ['../glads/00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = '00_pressure_seasonal_200_t190.png'
fig_00 = plot_pressure_maps_timeseries(fnames, figname, Qmin=10, Qmax=200, melt_forcing='SHMIP',
    t_ticklabels=t_ticklabels[:-1], t_xlabel=t_xlabel, t_ticks=t_ticks[:-1], t_lim=[1 + 3/12, 1 + 9/12],
    ff_ylim=[0,1.75], ff_yticks=[0, 0.5, 1, 1.5], tslice=190+365)


## Case 00a: Flat topo, scaled SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/00a_shimp_adj_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = '00_pressure_seasonal_scaled_melt_t190.png'
fig_00 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, melt_forcing='SHMIPadj',
    t_ticklabels=t_ticklabels[:-1], t_xlabel=t_xlabel, t_ticks=t_ticks[:-1], t_lim=[1 + 3/12, 1 + 9/12],
    tslice=190+365)


## Case 01: Flat topo, KAN_L forcing
KAN_tslice = 190+365
cases = [201, 202, 203, 204, 205]
pattern = '../glads/01_kan_l_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = '01_pressure_seasonal_t190.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, tslice=KAN_tslice, Qmin=1, Qmax=100, melt_forcing='KAN',
    t_ticklabels=t_ticklabels, t_xlabel=t_xlabel, t_ticks=t_ticks, t_lim=t_lim,
    ff_ylim=[0, 1.75], ff_yticks=[0, 0.5, 1, 1.5])


# ## Case 02: Trough topo, SHMIP forcing
# cases = [201, 202, 203, 204, 205]
# pattern = '../glads/02_shmip_forcing_synth_topo/RUN/output_%03d_seasonal.nc'
# fnames = [pattern % caseid for caseid in cases]
# figname = 'tmp.png'
# fig_02 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, melt_forcing='SHMIP',
#     t_ticklabels=t_ticklabels[:-1], t_xlabel=t_xlabel, t_ticks=t_ticks[:-1], t_lim=[1+3/12, 1+9/12],
#     ff_ylim=[0,1.75], ff_yticks=[0, 0.5, 1, 1.5], tslice=170+365)


## Case 02a: Trough topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
pattern = '../glads/02a_shmip_forcing_valley_topo/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = '02a_pressure_seasonal_t190.png'
fig_02 = plot_pressure_maps_timeseries(fnames, figname, Qmin=10, Qmax=200, melt_forcing='SHMIP',
    t_ticklabels=t_ticklabels[:-1], t_xlabel=t_xlabel, t_ticks=t_ticks[:-1], t_lim=[1+3/12, 1+9/12],
    ff_ylim=[0,1.75], ff_yticks=[0, 0.5, 1, 1.5], tslice=190+365)

"""
## Case 02: Valley topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
pattern = '../glads/02a_shmip_forcing_valley_topo/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = '02_pressure_seasonal_valley.png'
fig_02 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, melt_forcing='SHMIP',
    t_ticklabels=t_ticklabels[:-1], t_xlabel=t_xlabel, t_ticks=t_ticks[:-1], t_lim=[1+3/12, 1+9/12])


## Case 03: Trough topo, KAN_L forcing
cases = [101, 102, 103, 104, 105]
pattern = '../glads/03_kan_l_forcing_synth_topo/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = '03_pressure_seasonal.png'
fig_03 = plot_pressure_maps_timeseries(fnames, figname, tslice=KAN_tslice, Qmin=1, melt_forcing='KAN',
    t_ticklabels=t_ticklabels, t_xlabel=t_xlabel, t_ticks=t_ticks, t_lim=t_lim)
"""
plt.show()
