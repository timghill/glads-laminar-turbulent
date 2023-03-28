import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec


## Start the figure
figsize=(5, 4)
fig = plt.figure(figsize=figsize)


gs_kwargs=dict(wspace=0.05, hspace=0.2, 
        height_ratios=(2, 1), top=0.9)
gs = gridspec.GridSpec(2, 1, **gs_kwargs)

axs = np.array([fig.add_subplot(gs[i, 0]) for i in range(2)])

alphabet = ['a', 'b']
text_args = {'fontweight':'bold'}

## Read data
data = np.loadtxt('andrews2014_data.csv', skiprows=2, delimiter=',', usecols=(0, 4, 6), max_rows=5500)
tt_days = data[:, 0]
velocity = data[:, 1]
temp = data[:, 2]
tt_months = tt_days/365 * 12

print(tt_days)
print(tt_months)
print(velocity)

axs[0].plot(tt_months, velocity, linewidth=1, color='k')
axs[1].plot(tt_months, temp, linewidth=1, color='k')

xticks = np.array([4, 5, 6, 7, 8, 9, 10])
axs[0].set_xticks(xticks)
axs[0].set_xticklabels([])

axs[1].set_xticks(xticks)
axs[1].set_xticklabels(['May', '', 'July', '', 'Sep', '', 'Nov'])

axs[0].grid(linestyle=':', linewidth=0.5)
axs[1].grid(linestyle=':', linewidth=0.5)

plt.show()