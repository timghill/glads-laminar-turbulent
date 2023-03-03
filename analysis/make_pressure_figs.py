"""

Call plot_pressure_maps_timeseries.py to make pressure figures

"""

from matplotlib import pyplot as plt

from plot_pressure_maps_timeseries import plot_pressure_maps_timeseries

## Case 00: Flat topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/_00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = '00_pressure_seasonal.png'
fig_00 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1)

## Case 01: Flat topo, KAN_L forcing
KAN_tslice = 569
cases = [101, 102, 104, 104, 105]
pattern = '../glads/01_kan_l_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = '01_pressure_seasonal.png'
fig_01 = plot_pressure_maps_timeseries(fnames, figname, tslice=KAN_tslice, Qmin=1)

## Case 03: Trough topo, KAN_L forcing
cases = [101, 102, 104, 104, 105]
pattern = '../glads/03_kan_l_forcing_synth_topo/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = '03_pressure_seasonal.png'
fig_03 = plot_pressure_maps_timeseries(fnames, figname, tslice=KAN_tslice, Qmin=1)

plt.show()
