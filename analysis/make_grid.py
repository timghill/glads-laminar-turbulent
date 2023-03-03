"""

Call plot_pressure_maps_timeseries.py to make pressure figures

"""

from matplotlib import pyplot as plt

from plot_pressure_grid import plot_pressure_grid

## Case 00: Flat topo, SHMIP forcing
cases = [[1, 2, 3, 4, 5],
         [101, 102, 104, 104, 105],
         [201, 202, 203, 204, 205],
         [101, 102, 104, 104, 105]
]

patterns = ['../glads/_00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc',
            '../glads/01_kan_l_forcing_shmip_topo/RUN/output_%03d_seasonal.nc',
            '../glads/02_shmip_forcing_synth_topo/RUN/output_%03d_seasonal.nc',
            '../glads/03_kan_l_forcing_synth_topo/RUN/output_%03d_seasonal.nc',
]

fnames = [[patterns[j] % cases[j][i] for i in range(5)] for j in range(4)]
print(fnames)
figname = 'pressure_grid.png'
fig_00 = plot_pressure_grid(fnames, figname)

plt.show()
