import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt

proj_cases = [1, 3]
proj_fnames = ['../RUN/output_%03d_seasonal.nc' % caseid for caseid in proj_cases]

home_cases = [301, 302]
home_fnames = ['/home/tghill/subglacial-emulator/glads/seasonal_KAN/RUN/output_%03d_seasonal.nc' % caseid for caseid in home_cases]

figname = 'S_comparison.png'

#cases = [101, 102, 103]
#fname_pattern = '../RUN/output_%03d_seasonal.nc'
#figname = 'project_channels_10x.png'

#cases = [401, 402]
#fname_pattern = '/home/tghill/scratch/subglacial-emulator/glads/seasonal_KAN/RUN/output_%03d_seasonal.nc'
#figname = 'scratch_channels.png'

colors= ['r', 'm', 'b', 'c']

fig, ax = plt.subplots()
for i in range(len(proj_fnames)):
    fname = proj_fnames[i]
    print(fname)
    out = nc.Dataset(fname)
    S = out['S_channel'][:].data.T
    time = out['time'][:].data
    time = (time - time[0])/86400/365
    labels = ['Proj Turb', 'Proj Lam']
    colors = ['sandybrown', 'firebrick']
    ax.plot(time, np.sum(S, axis=0), color=colors[i], label=labels[i])

for i in range(len(home_fnames)):
    fname = home_fnames[i]
    print(fname)
    out = nc.Dataset(fname)
    S = out['S_channel'][:].data.T
    time = out['time'][:].data
    time = (time - time[0])/86400/365
    labels=['Home Turb', 'Home Lam']
    colors = ['cyan', 'slateblue']
    ax.plot(time, np.sum(S, axis=0), color=colors[i], label=labels[i])

ax.legend()
fig.savefig(figname, dpi=600)
plt.show()

