import sys
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as apc
import astropy.units as apu

sys.path.append('/Users/amsys/py4py')
from py4py import load_grid, plot_dat_many


x, z = load_grid('qso_100.dom0.')
titles = [r'$\alpha=0.6$, R$_v$=10$^{19}$',
          r'$\alpha=1.0$, R$_v$=10$^{18}$']
xlims = [(16.6, 19.0), (16.6, 19.0)]
zlims = [(14.9, 18.5), (14.9, 18.5)]
#
# table_ion_h1_normal = Table.read('qso_100.ionH1.dat', format='ascii')
# table_ion_h1_tweak = Table.read('qso_100_tweak.ionH1.dat', format='ascii')
# table_ion_c4_normal = Table.read('qso_100.ionC4.dat', format='ascii')
# table_ion_c4_tweak = Table.read('qso_100_tweak.ionC4.dat', format='ascii')
#
# table_line_ha_normal = Table.read('qso_100.lineH1.3-2.dat', format='ascii')
# table_line_ha_tweak = Table.read('qso_100_tweak.lineH1.3-2.dat', format='ascii')
# table_line_c4_normal = Table.read('qso_100.lineC4.dat', format='ascii')
# table_line_c4_tweak = Table.read('qso_100_tweak.lineC4.dat', format='ascii')
#
# table_density_normal = Table.read('qso_100.rho.dat', format='ascii')
# table_density_tweak = Table.read('qso_100_tweak.rho.dat', format='ascii')
#
# fig_ion_h1 = plot_dat_many([table_ion_h1_normal, table_ion_h1_tweak], [x, x], [z, z],
#               xlims=xlims, zlims=zlims, titles=titles,
#               title='H-I ion fraction', label='Log ion fraction',
#               shared_y=True, volume=False, shared_cbar=True)
# fig_ion_h1.savefig('ion_h1.eps')
#
# fig_ion_c4 = plot_dat_many([table_ion_c4_normal, table_ion_c4_tweak], [x, x], [z, z],
#               xlims=xlims, zlims=zlims, titles=titles,
#               title='C-IV ion fraction', label='Log ion fraction',
#               shared_y=True, volume=False, shared_cbar=True)
# fig_ion_c4.savefig('ion_c4.eps')
#
# fig_line_c4 = plot_dat_many([table_line_c4_normal, table_line_c4_tweak], [x, x], [z, z],
#               xlims=xlims, zlims=zlims, titles=titles,
#               title='rC$_{IV}$ luminosity', label='Log luminosity (erg cm$^{-3}$ s$^{-1}$)',
#               shared_y=True, volume=True, shared_cbar=True)
# fig_line_c4.savefig('line_c4.eps')
#
# fig_line_ha = plot_dat_many([table_line_ha_normal, table_line_ha_tweak], [x, x], [z, z],
#               xlims=xlims, zlims=zlims, titles=titles,
#               title=r'H$\alpha$ luminosity', label='Log luminosity (erg cm$^{-3}$ s$^{-1}$)',
#               shared_y=True, volume=True, shared_cbar=True)
# fig_line_ha.savefig('line_ha.eps')
#
# fig_rho = plot_dat_many([table_density_normal, table_density_tweak], [x, x], [z, z],
#               xlims=xlims, zlims=zlims, titles=titles,
#               title=r'Wind density', label=r'Density $\rho$ (g cm$^{-3}$)',
#               shared_y=True, volume=False, shared_cbar=True)
# fig_rho.savefig('rho.eps')
#
# fig_rho_lin = plot_dat_many([table_density_normal, table_density_tweak], [x, x], [z, z],
#               xlims=xlims, zlims=zlims, titles=titles,
#               title=r'Wind density', label=r'Density $\rho$ (g cm$^{-3}$)',
#               shared_y=True, volume=False, shared_cbar=True, log=False)
# fig_rho_lin.savefig('rho_linear.eps')


k_low = np.zeros(100)
k_one = np.zeros(100)
v_low_inner = np.zeros(100)
v_one_inner = np.zeros(100)
v_low_outer = np.zeros(100)
v_one_outer = np.zeros(100)

cm = apu.m / 100

t = Table(
    [k_low, k_one, v_low_inner,  v_low_outer, v_one_inner, v_one_outer],
    names=['k_low', 'k_one', 'v_low_inner', 'v_low_outer', 'v_one_inner', 'v_one_outer']
)
t['v_low_inner'].unit = cm / apu.second
t['v_low_outer'].unit = cm / apu.second
t['v_one_inner'].unit = cm / apu.second
t['v_one_outer'].unit = cm / apu.second

v_zero = 6e5 * cm / apu.second
v_escape_inner = np.sqrt(2. * apc.G * apc.M_sun * 1e9 / (8.85667e14 * cm))
v_escape_outer = np.sqrt(2. * apc.G * apc.M_sun * 1e9 / (5e17 * cm))
v_inf_inner = 1 * v_escape_inner
v_inf_outer = 1 * v_escape_outer

r = np.logspace(15, 19, num=100)
for i, r_i in enumerate(r):
    t['k_low'][i] = np.power(r_i/1e19, 0.6) / (1 + np.power(r_i/1e19, 0.6))
    t['k_one'][i] = (r_i/1e18) / (1 + r_i/1e18)
    t['v_low_inner'].quantity[i] = v_zero + (v_inf_inner - v_zero) * t['k_low'][i]
    t['v_one_inner'].quantity[i] = v_zero + (v_inf_inner - v_zero) * t['k_one'][i]
    t['v_low_outer'].quantity[i] = v_zero + (v_inf_outer - v_zero) * t['k_low'][i]
    t['v_one_outer'].quantity[i] = v_zero + (v_inf_outer - v_zero) * t['k_one'][i]

fig, ax = plt.subplots(1)
ax.set_xlabel('Streamline position (cm)')
ax.set_ylabel('Velocity (cm/s)')
ax.plot(np.log10(r), t['v_low_inner'], '-r', label=titles[0])
ax.plot(np.log10(r), t['v_one_inner'], '-b', label=titles[1])
ax.legend()
fig.savefig('v_factor_inner.eps')

fig, ax = plt.subplots(1)
ax.set_xlabel('Streamline position (cm)')
ax.set_ylabel('Velocity (cm/s)')
ax.fill_between(np.log10(r), t['v_low_inner'], t['v_low_outer'], facecolor=(1,0,0,0.25), edgecolor=(1,0,0), label=titles[0])
ax.fill_between(np.log10(r), t['v_one_inner'], t['v_one_outer'], facecolor=(0,0,1,0.25), edgecolor=(0,0,1), label=titles[1])
ax.legend(loc='upper left')
fig.savefig('v_factor_both.png')
