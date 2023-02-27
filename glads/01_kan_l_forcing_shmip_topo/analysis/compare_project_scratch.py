import numpy as np
from matplotlib import pyplot as plt
import netCDF4 as nc

project = nc.Dataset('../RUN/output_001_seasonal.nc');
scratch = nc.Dataset('/home/tghill/subglacial-emulator/glads/seasonal_KAN/RUN/output_301_seasonal.nc')

# fig, (ax1, ax2) = plt.subplots()

"""
    variables(dimensions): float64 rho_w(), float64 rho_i(), float64 g_grav(), float64 L_fusion(), float64 c_w(), float64 c_t_c(), float64 hour(), float64 day(), float64 year(), float64 n_glen(), float64 creep_const_s(), float64 creep_const_s_soft(), float64 creep_const_s_trans(), float64 creep_const_c(), float64 l_bed(), float64 l_c(), float64 h_bed(), float64 alpha_c(), float64 beta_c(), float64 cond_c(), float64 nu(), float64 omega(), float64 cond_s(cond), float64 alpha_s(cond), float64 beta_s(cond)
"""

print(project['para/cond_s'][:].data)
print(scratch['para/cond_s'][:].data)

print(project['para/l_c'][:].data)
print(scratch['para/l_c'][:].data)

print(project['para/cond_c'][:].data)
print(scratch['para/cond_c'][:].data)

print(project['para/alpha_c'][:].data)
print(scratch['para/alpha_c'][:].data)

proj_bed = project['bed'][:].data
scra_bed = scratch['bed'][:].data

print(proj_bed - scra_bed)

proj_thick = project['thick'][:].data
scra_thick = scratch['thick'][:].data

print(proj_thick - scra_thick)
