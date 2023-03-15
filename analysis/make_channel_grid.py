"""

Call plot_channel_volume.py to make pressure figures

"""

import numpy as np

from matplotlib import pyplot as plt

from plot_pressure_grid import plot_pressure_grid

## Global

t_lim = [1 + 4/12, 1 + 10/12]
t_ticks = 1 + np.arange(4, 11)/12
t_ticklabels = ['4', '', '6', '', '8', '', '10']
ff_ylim = [0, 1.75]

cases = [[201, 202, 203, 204, 205],
         [101, 102, 103, 104, 105],
         [301, 302, 303, 304, 305],
         [101, 102, 103, 104, 105]
]

patterns = ['../glads/00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc',
            '../glads/01_kan_l_forcing_shmip_topo/RUN/output_%03d_seasonal.nc',
            '../glads/02_shmip_forcing_synth_topo/RUN/output_%03d_seasonal.nc',
            '../glads/03_kan_l_forcing_synth_topo/RUN/output_%03d_seasonal.nc',
]

fnames = [[patterns[j] % cases[j][i] for i in range(5)] for j in range(4)]
figname = 'pressure_grid_new.png'
fig_00 = plot_channel_volume(fnames, figname,
    tlim=t_lim, t_ticks=t_ticks, t_ticklabels=t_ticklabels,
    xlabel='Month')

plt.show()
