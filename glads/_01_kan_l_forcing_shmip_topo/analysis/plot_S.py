import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt

cases = [101, 102, 103]
fname_pattern = '../RUN/output_%03d_seasonal.nc'
figname = 'project_channels_10x.png'

# cases = [301, 302, 303]
# fname_pattern = '/home/tghill/scratch/subglacial-emulator/glads/seasonal_KAN/RUN/output_%03d_seasonal.nc'
# figname = 'scratch_channels.png'

fig, ax = plt.subplots()
for i in range(len(cases)):
    fname = fname_pattern % cases[i]
    out = nc.Dataset(fname)
    S = out['S_channel'][:].data.T
    time = out['time'][:].data
    time = (time - time[0])/86400/365

    ax.plot(time, np.sum(S, axis=0))

fig.savefig(figname, dpi=600)
plt.show()

