"""

Compare summer and winter water pressure for turbulent, transition,
and laminar models

"""

import numpy as np
import netCDF4 as nc

from matplotlib import pyplot as plt

# Global constants
rhow = 1000
rhoi = 910
g = 9.81

x_band = [70e3]
band_width = 60e3

x_bands = np.array([15e3, 30e3, 70e3])
band_widths = 5e3
winter_tindex = [350, 365]
labels = ['Turbulent 5/4', 'Turbulent 3/2',
          'Laminar',
          'Transition 5/4', 'Transition 3/2']
labels = [l.ljust(14) for l in labels]

def quantify_floatation(fnames, metric=np.nanmean, tslices=None,
    x_bands=None, band_width=None, strfmt='%.4f', labels=None):
    """Docstring..."""
    if x_bands is None:
        ff_means = np.zeros((len(fnames), 1))
    else:
        ff_means = np.zeros((len(fnames), len(x_bands)))

    if tslices is None:
        tslices = [0, -1]

    if labels is None:
        labels = [fname.split('/')[-1] for fname in fnames]

    ff_strfmt = []
    for i in range(len(fnames)):
        out = nc.Dataset(fnames[i])
        phi = out['phi'][tslices[0]:tslices[1], :].data.T
        N = out['N'][tslices[0]:tslices[1], :].data.T
        
        bed = np.vstack(out['bed'][:].data)
        phi_elevation = rhow*g*bed

        pw = phi - phi_elevation
        ff = pw/(N + pw)
        
        if x_bands is None:
            ff_mean = metric(ff)
            ff_means[i, 0] = ff_mean
            ff_strfmt.append(labels[i] + ':\t' + (strfmt % ff_mean))
        else:
            for k in range(len(x_bands)):
                xmin = x_bands[k] - band_width/2
                xmax = x_bands[k] + band_width/2
                nodex = out['nodes'][0, :].data
                mask = np.logical_and(nodex>=xmin, nodex<=xmax)
                ff_mean = metric(ff[mask, :])
                ff_means[i, k] = ff_mean
            ff_strfmt.append(labels[i] + ':\t' + '\t'.join(
                    [strfmt % f for f in ff_means[i, :]]))
        out.close()

    ff_strfmt = '\n'.join(ff_strfmt)
    return ff_means, ff_strfmt


print('-------------------------------------------------')
print('    STEADY (case 00)')
print('-------------------------------------------------')

steady_cases = [1, 2, 3, 4, 5]
steady_dir = '../glads/00_shmip_forcing_shmip_topo/RUN/output_%03d_steady.nc'
fnames = [steady_dir % caseid for caseid in steady_cases]

print('\nMean water pressure')
data, strfmt = quantify_floatation(fnames, x_bands=x_bands, band_width=band_widths,
    tslices=[-2, -1], labels=labels)
print(strfmt)

print('\nPeak water pressure')
data, strfmt = quantify_floatation(fnames, x_bands=x_bands, band_width=band_widths,
    metric=np.nanmax, labels=labels, tslices=[-2, -1])
print(strfmt)


print('\n\n')
print('-------------------------------------------------')
print('    SEASONAL (case 00)')
print('-------------------------------------------------')

seasonal_cases = [1, 2, 3, 4, 5]
seasonal_dir = '../glads/00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'
fnames = [seasonal_dir % caseid for caseid in seasonal_cases]

print('\nMean winter water pressure')
data, strfmt = quantify_floatation(fnames, x_bands=x_bands, band_width=band_widths,
    tslices=winter_tindex, labels=labels)
print(strfmt)

print('\nPeak summer water pressure')
data, strfmt = quantify_floatation(fnames, x_bands=x_bands, band_width=band_widths,
    metric=np.nanmax, labels=labels)
print(strfmt)

print('\n\n')
print('-------------------------------------------------')
print('    SEASONAL (case 01)')
print('-------------------------------------------------')

seasonal_cases = [1, 1, 1, 1, 1]
seasonal_dir = '../glads/01_kan_l_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'
fnames = [seasonal_dir % caseid for caseid in seasonal_cases]

print('\nMean winter water pressure')
data, strfmt = quantify_floatation(fnames, x_bands=x_bands, band_width=band_widths,
    tslices=winter_tindex, labels=labels)
print(strfmt)

print('\nPeak summer water pressure')
data, strfmt = quantify_floatation(fnames, x_bands=x_bands, band_width=band_widths,
    metric=np.nanmax, labels=labels)
print(strfmt)


"""
xb_print = list(x_bands/1e3)
print('Average winter floatation fraction for bands: ', xb_print)
for i in range(len(seasonal_cases)):
    out = nc.Dataset(seasonal_dir % seasonal_cases[i])
    phi = out['phi'][winter_index[0]:winter_index[1], :].data.T
    N = out['N'][winter_index[0]:winter_index[1], :].data.T
    bed = np.vstack(out['bed'][:].data)
    phi_elevation = rhow*g*bed
    pw = phi - phi_elevation
    ff = pw/(N + pw)
    
    ff_means = len(seasonal_cases)*['']
    for k in range(len(x_bands)):
        xmin = x_bands[k] - band_width/2
        xmax = x_bands[k] + band_width/2
        nodex = out['nodes'][0, :].data
        mask = np.logical_and(nodex>=xmin, nodex<=xmax)
        ff_mean = np.nanmean(ff[mask, :])
        ff_means[k] = '%.4f' % ff_mean
    print(labels[i]+':\t' + '\t'.join(ff_means))

    out.close()


print('\nPeak summer floatation fraction for bands: ', xb_print)
for i in range(len(seasonal_cases)):
    out = nc.Dataset(seasonal_dir % seasonal_cases[i])
    phi = out['phi'][:].data.T
    N = out['N'][:].data.T
    bed = np.vstack(out['bed'][:].data)
    phi_elevation = rhow*g*bed
    pw = phi - phi_elevation
    ff = pw/(N + pw)
    
    ff_peaks = len(seasonal_cases)*['']
    for k in range(len(x_bands)):
        xmin = x_bands[k] - band_width/2
        xmax = x_bands[k] + band_width/2
        nodex = out['nodes'][0, :].data
        mask = np.logical_and(nodex>=xmin, nodex<=xmax)
        ff_peak = np.nanmax(ff[mask, :])
        ff_peaks[k] = '%.4f' % ff_peak
    print(labels[i]+':\t' + '\t'.join(ff_peaks))
    out.close()



print('\n\n')
print('-------------------------------------------------')
print('    SEASONAL (case 01)\n')

# Define winter in terms of the range of time indices
winter_index = [350, 365]

# Range of domain to use
x_bands = np.array([15e3, 30e3, 70e3])
band_width = 5e3

seasonal_cases = [1, 1, 1, 1, 1]

labels = ['Turbulent 5/4', 'Turbulent 3/2',
          'Laminar', 'Transition 5/4', 'Transition 3/2']
labels = [l.ljust(14) for l in labels]
seasonal_dir = '../glads/01_kan_l_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'

xb_print = list(x_bands/1e3)
print('Average winter floatation fraction for bands: ', xb_print)
for i in range(len(seasonal_cases)):
    out = nc.Dataset(seasonal_dir % seasonal_cases[i])
    phi = out['phi'][winter_index[0]:winter_index[1], :].data.T
    N = out['N'][winter_index[0]:winter_index[1], :].data.T
    bed = np.vstack(out['bed'][:].data)
    phi_elevation = rhow*g*bed
    pw = phi - phi_elevation
    ff = pw/(N + pw)
    
    ff_means = len(seasonal_cases)*['']
    for k in range(len(x_bands)):
        xmin = x_bands[k] - band_width/2
        xmax = x_bands[k] + band_width/2
        nodex = out['nodes'][0, :].data
        mask = np.logical_and(nodex>=xmin, nodex<=xmax)
        ff_mean = np.nanmean(ff[mask, :])
        ff_means[k] = '%.4f' % ff_mean
    print(labels[i]+':\t' + '\t'.join(ff_means))
    out.close()


print('\nPeak summer floatation fraction for bands: ', xb_print)
for i in range(len(seasonal_cases)):
    out = nc.Dataset(seasonal_dir % seasonal_cases[i])
    phi = out['phi'][:].data.T
    N = out['N'][:].data.T
    bed = np.vstack(out['bed'][:].data)
    phi_elevation = rhow*g*bed
    pw = phi - phi_elevation
    ff = pw/(N + pw)

    ff_peaks = len(seasonal_cases)*['']
    for k in range(len(x_bands)):
        xmin = x_bands[k] - band_width/2
        xmax = x_bands[k] + band_width/2
        nodex = out['nodes'][0, :].data
        mask = np.logical_and(nodex>=xmin, nodex<=xmax)
        ff_peak = np.nanmax(ff[mask, :])
        ff_peaks[k] = '%.4f' % ff_peak
    print(labels[i]+':\t' + '\t'.join(ff_peaks))
    out.close()
"""







