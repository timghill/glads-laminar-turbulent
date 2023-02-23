"""
Temporary fix for lapse rate for KAN_L temperature record.
"""

import numpy as np
from matplotlib import pyplot as plt

temps = np.loadtxt('KAN_L_2014_temp_clipped.txt', delimiter=',')
time = temps[:, 0]
T = temps[:, 1]

z_stn = 665
lr_old = -0.0075
lr_new = -0.005

T_stn = T + lr_old*z_stn

T_sl = T_stn + lr_new*(0 - z_stn)

fig, ax = plt.subplots()
ax.plot(time, T, label=r'SL Temp ($\Gamma = -0.0075^{\circ}\rm{C}~\rm{m}^{-1}$)')
ax.plot(time, T_sl, label=r'SL Temp ($\Gamma = -0.005^{\circ}\rm{C}~\rm{m}^{-1}$)')
ax.plot(time, T_stn, label=r'KAN_L Temp')
ax.legend()

temps_updated = np.zeros(temps.shape)
temps_updated[:, 0] = time
temps_updated[:, 1] = T_sl

np.savetxt('KAN_L_2014_temp_clipped_lr.txt', temps_updated, fmt=('%d, %.3e'))

plt.show()
