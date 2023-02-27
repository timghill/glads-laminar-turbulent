"""

Investigate mechanisms behind differences for turbulent/laminar models

"""

import numpy as np
import netCDF4 as nc

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from scipy import interpolate

# File structure and cases
fname_pattern = '../RUN/output_%03d_seasonal.nc'
cases = [911, 911, 912, 913, 913]
figname = 'mechanisms.png'
labels = ['Turbulent 5/4', 'Turbulent 3/2', 'Laminar', 'Transition 5/4', 'Transition 3/2']
textx = 0.1
texty = 0.85
textfmt = {'fontweight':'bold', 'fontsize':10}

tband_min = 20e3
tband_max = 60e3


# Get filenames and num cases
n_cases = len(cases)

# Constants
k_turb = 0.00919703065356383
rhow = 1000
g = 9.81

# Helper function
def width_average(xy, z, dx=1e3):
    """Width-average field z according to coordinates xy.

    Arguments:
    ---------
    xy : (N, 2) array in meters
    z : (N,) array
    dx : Flowline increment in meters
    """
    xedge = np.arange(0, 100e3+dx, dx)
    xmid = 0.5*(xedge[1:] + xedge[:-1])

    xvec = np.reshape(xy[:, 0], (xy.shape[0], 1))
    gvec = np.array([xmid])

    tf_mask = np.logical_and(xvec>=(gvec-dx/2), xvec<=(gvec+dx/2))

    z_avg = np.zeros(xmid.shape)
    for i in range(len(xmid)):
        z_avg[i] = np.mean(z[tf_mask[:, i]])
    
    return xmid, z_avg

# Set figure
fig = plt.figure(figsize=(6, 6))
gs = GridSpec(5, 2, hspace=0.05, top=0.9, bottom=0.075, right=0.95, left=0.15)
alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
axs = np.array([[fig.add_subplot(gs[ii, jj]) for jj in range(2)] for ii in range(5)])

dx = 1
xgrid = np.arange(0, 100, dx)

w = 25e3    # Domain width (km)

colors = np.array([[0.420, 0.510, 0.620, 1],
                   [0.579, 0.677, 0.781, 1],
                   [0.500, 0.500, 0.500, 1],
                   [0.859, 0.683, 0.275, 1],
                   [0.929, 0.835, 0.408, 1]])

for ii in range(n_cases):
    fname = fname_pattern % cases[ii]
    print(fname)

    out = nc.Dataset(fname, 'r')

    nodes = out['nodes'][:].data.T
    elements = out['elements'][:].data.T
    connect_edge = out['connect_edge'][:2, :].data.T.astype(int) - 1

    band_mask = np.logical_and(elements[:, 0]>=tband_min, elements[:, 1]<=tband_max)

    n_nodes = nodes.shape[0]
    n_elements = elements.shape[0]
    n_edges = connect_edge.shape[0]

    time = out['time'][:].data.T
    n_time = len(time)

    ## Net sheet and channel flux
    qsx = out['qs'][:, 0, :].data.T
    qsy = out['qs'][:, 1, :].data.T
    qs_xy = out['qs'][:].data.T
    qs = np.sqrt(qs_xy[:, 0]**2 + qs_xy[:, 1]**2)

    Q = out['Q'][:].data.T

    Qtot = np.zeros((xgrid.shape[0], n_time))
    qtot = np.zeros((xgrid.shape[0], n_time))
    for jj in range(len(xgrid)-1):
        x1 = xgrid[jj]
        x2 = x1 + dx

        # Hack -- this is not quite right
        el_mask = np.logical_and(elements[:, 0]/1e3>=x1, elements[:, 0]/1e3<x2)
        qsum = w*np.nanmean(qs[el_mask, :], axis=0)

        n1 = nodes[connect_edge[:, 0], :]
        n2 = nodes[connect_edge[:, 1], :]
        edge_mask = np.logical_and(n1[:, 0]/1e3>=x1, n2[:, 0]/1e3<x2)
        Qsum = np.sum(np.abs(Q[edge_mask, :]), axis=0)

        qtot[jj, :] = qsum
        Qtot[jj, :] = Qsum

    ax1 = axs[0, 0]
    ax2 = axs[0, 1]

    tt = time/86400/365 - 100
    ax1.plot(xgrid, qtot[:, 183 + 365], color=colors[ii])
    # ax1.plot(xgrid/1e3, Qtot[:, 183 + 365], color=colors[ii], linestyle='-.')
    ax1.grid()
    ax1.set_ylim([0, 1000])
    ax1.set_ylabel(r'q (m$^3$ s$^{-1}$)')
    ax1.set_xlim([0, 100])
    ax1.legend(labels=labels, bbox_to_anchor=[0., 1, 2.2, 0.2], ncol=3,
        loc='lower center', frameon=False, mode='expand')
    ax1.text(textx, texty, 'a', transform=ax1.transAxes, **textfmt)


    ax2.plot(tt, np.mean(qtot[np.logical_and(xgrid>=tband_min/1e3, xgrid<=tband_max/1e3), :], axis=0), color=colors[ii])
    # ax2.plot(tt, np.mean(Qtot, axis=0), color=colors[ii], linestyle='-.')
    ax2.grid()
    ax2.set_ylim([0, 400])
    ax2.set_xlim([1, 2])
    ax2.set_xticks([1, 1.25, 1.5, 1.75, 2])
    ax2.text(textx, texty, 'b', transform=ax2.transAxes, **textfmt)
    ## Sheet transmissivity
    # Sheet flux qs:        elements
    # Sheet thickness h:    nodes
    # Sheet pot grad:       elements
    # --> interpolate thickness onto elements
    hs = out['h_sheet'][:].data.T
    interpolator = interpolate.LinearNDInterpolator(nodes, hs)
    hs_element = interpolator(elements)
    k = float(out['para/cond_s'][:].data)

    pp = 1
    if ii==0 or ii==1:
        gradphi = (qs/k/hs_element**(5/4))**2
        
        # T = k*hs_element**(5/4)*gradphi**(-1/2)
        # T = rhow*g*k*hs_element**(5/4)*gradphi**(-1/2)

        # hpow = k**pp * hs_element**(5/4)

        # phipow = k**(1-pp)*gradphi**(1/2)
    elif ii==2:
        gradphi = (qs/k/hs_element**(3))**1
        # T = k*hs_element**3
        # T = k*hs_element**3
        # hpow = k**pp * hs_element**(3)
        # phipow = k**(1-pp) * gradphi**(1)
    else:
        omega = float(out['para/omega'][:].data)
        nu = float(out['para/nu'][:].data)
        h_bed = float(out['para/h_bed'][:].data)
        gradphi = ((qs + omega/nu * (hs_element/h_bed)**(1/2) * qs**2)/k/hs_element**(3))**1

        # sq_arg = 1 + 4*omega/nu*(hs_element/h_bed)**(0.5)*k*hs_element**3*gradphi
        # T = nu*np.sqrt(h_bed/hs_element)/omega * (-1 + np.sqrt(sq_arg)) / gradphi
        # hpow = k**pp * hs_element**(3)
        # phipow = k**(1-pp) * gradphi**(1)

        # T = k*hs_element**3

    T = rhow*g*(qs/gradphi)
    
    x_mean, T_mean = width_average(elements, T[:, 183+365])
    ax1 = axs[1, 0]
    ax1.plot(x_mean/1e3, T_mean, color=colors[ii])
    # ax1.scatter(elements[:, 0], T[:, 183+365], 1, color=colors[ii])
    ax1.set_ylim([0, 5])
    # ax1.grid()
    # ax1.pcolormesh(T, vmin=0, vmax=0.1)
    ax1.set_ylabel(r'T (m$^2$ s$^{-1}$)')
    ax1.grid()
    ax1.set_xlim([0, 100])
    ax1.text(textx, texty, 'c', transform=ax1.transAxes, **textfmt)

    T[elements[:, 0]<10e3, :] = np.nan
    ax2 = axs[1, 1]
    # T_tt = np.nanmean(T, axis=0)
    # if ii>0:
    ax2.plot(tt, np.nanmean(T[band_mask, :], axis=0), color=colors[ii])
    ax2.grid()
    ax2.set_xlim([1, 2])
    ax2.set_xticks([1, 1.25, 1.5, 1.75, 2])
    ax2.text(textx, texty, 'd', transform=ax2.transAxes, **textfmt)
    # ax2.plot(tt, np.mean(hs_element, axis=0), color=colors[ii])

    ## h to appropriate power
    ax1 = axs[2, 0]
    # ax1.scatter(elements[:, 0], hpow[:, 183+365], 1, color=colors[ii])
    x_mean, hs_mean = width_average(elements, hs_element[:, 183+365])
    ax1.plot(x_mean/1e3, hs_mean, color=colors[ii])
    # ax1.scatter(elements[:, 0], hs_element[:, 183+365], 1, color=colors[ii])
    ax1.set_ylabel(r'$h$ (m)')
    ax1.set_ylim([0, 0.22])
    ax1.grid()
    ax1.set_xlim([0, 100])
    ax1.text(textx, texty, 'e', transform=ax1.transAxes, **textfmt)
    
    # hpow[elements[:, 0]<1e3, :] = np.nan
    ax2 = axs[2, 1]
    ax2.plot(tt, np.nanmean(hs_element[band_mask, :], axis=0), color=colors[ii])
    ax2.grid()
    ax2.set_xlim([1, 2])
    ax2.set_xticks([1, 1.25, 1.5, 1.75, 2])
    ax2.text(textx, texty, 'f', transform=ax2.transAxes, **textfmt)

    ## gradphi to power
    ax1 = axs[3, 0]
    x_mean, gradphi_mean = width_average(elements, gradphi[:, 183+365], dx=2e3)
    ax1.plot(x_mean/1e3, gradphi_mean, color=colors[ii])
    # ax1.scatter(elements[:, 0], gradphi[:, 183+365], 1, color=colors[ii])
    ax1.set_ylabel(r'$|\nabla \phi|$ (Pa m$^{-1}$)')
    ax1.grid()
    ax1.set_xlim([0, 100])
    ax1.text(textx, texty, 'g', transform=ax1.transAxes, **textfmt)
    # ax1.set_ylim([0, 600])

    ax2 = axs[3, 1]
    ax2.plot(tt, np.nanmean(gradphi[band_mask, :], axis=0), color=colors[ii])
    ax2.grid()
    ax2.set_xlim([1, 2])
    ax2.set_xticks([1, 1.25, 1.5, 1.75, 2])
    ax2.text(textx, texty, 'h', transform=ax2.transAxes, **textfmt)

    ## Effective turbulent conductivity
    k_eff = qs/hs_element**(5/4)/gradphi**(1/2)

    ax1 = axs[4, 0]
    x_mean, k_eff_mean = width_average(elements, k_eff[:, 183+365])
    ax1.plot(x_mean/1e3, k_eff_mean/k_turb, color=colors[ii])
    # ax1.scatter(elements[:, 0], k_eff[:, 183+365], 1, color=colors[ii])
    ax1.grid()
    ax1.set_ylabel(r'$k_{\rm{eff}}/k_{\rm{turb}}$')
    ax1.set_xlim([0, 100])
    ax1.text(textx, texty, 'i', transform=ax1.transAxes, **textfmt)

    ax2 = axs[4, 1]
    ax2.plot(tt, np.nanmean(k_eff[band_mask, :], axis=0)/k_turb, color=colors[ii])
    ax2.grid()
    ax2.set_xlim([1, 2])
    ax2.set_xticks([1, 1.25, 1.5, 1.75, 2])
    ax2.text(textx, texty, 'j', transform=ax2.transAxes, **textfmt)

ax1.set_xlabel('x (km)')
ax2.set_xlabel('Year')

for ax in fig.get_axes():
    ss = ax.get_subplotspec()
    ax.spines.top.set_visible(False)
    ax.spines.bottom.set_visible(True)
    ax.spines.left.set_visible(True)
    ax.spines.right.set_visible(False)

    ax.tick_params(axis='x', which='both', labelbottom=ss.is_last_row())
    ax.tick_params(axis='y', which='both', labeltop=ss.is_first_col())


fig.savefig(figname, dpi=600)
plt.show()

    
