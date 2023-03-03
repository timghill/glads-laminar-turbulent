"""

Call plot_mechanisms.py

"""

from plot_mechanisms import plot_mechanisms
import defaults

## Case 00: Flat topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/_00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = '00_mechanisms.png'
models = defaults.labels
fig_00 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.0137799761144865)

# ## Case 01: Flat topo, KAN_L forcing
# cases = [101, 102, 104, 104, 105]
# pattern = '../glads/01_kan_l_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'
# fnames = [pattern % caseid for caseid in cases]
# figname = '01_pressure_seasonal.png'
# fig_01 = plot_pressure_maps_timeseries(fnames, figname)