import seaborn as sns
from scipy import constants

hex_codes = sns.color_palette().as_hex()
hex_codes2 = sns.color_palette('husl', 9).as_hex()

from matplotlib.colors import ListedColormap

#flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
my_cmap = ListedColormap(hex_codes)
my_cmap2 = ListedColormap(hex_codes2)
my_cmap3 = ListedColormap(["#9b59b6", "#34495e", "#3498db", "#95a5a6", "#e74c3c"])

### example solutions

from exact_coupled_transmission_line_eqn_solver import *
from analytic_coupled_lines_calculator import calc_self_and_coupling_capacitance, calc_self_and_coupling_inductance

import cap_util as cap

w = 5e-6
s = 7.5e-6
eps_r = 11.7

l_Rn = 974e-6
l_Rf = 1617e-6
l_Gn = 759e-6
l_Gf = 1659e-6
l_c = 318e-6

d_val = 5.5*1e-6

cs, cm = calc_self_and_coupling_capacitance(d_val, w, s, eps_r)
ls, lm = calc_self_and_coupling_inductance(d_val, w, s)

Cm = cm[0]
Lm = lm[0]

large_dist = 200*1e-6
cs_large_dist, cm_large_dist = calc_self_and_coupling_capacitance(large_dist, w, s, eps_r)
ls_large_dist, lm_large_dist = calc_self_and_coupling_inductance(large_dist, w, s)

Z0 = np.sqrt(ls_large_dist/cs_large_dist)
phase_vel = np.sqrt(1/(ls_large_dist * cs_large_dist))

Z0 = Z0[0]
phase_vel = phase_vel[0]

print('Cm:', Cm)
print('Lm:', Lm)

print('Z0:', Z0)
print('phase_vel:', phase_vel)

L_readout = l_Gn + l_c + l_Gf
L_purcell = l_Rn + l_c + l_Rf

J_val = J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=phase_vel, Z0=Z0)
print('J_val_dist (MHz):', J_val/(2*np.pi*1e6))

omegas = np.arange(7.5115, 11.8115, 0.005) * 1e9 * 2*np.pi
test_Z_transfer_exact = Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

plt.figure(figsize=(8, 3.85))
test_Z_equiv_LE_circuit_symbolic = Z_transfer_equivalent_LE_circuit_from_symbolic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_exact), color = my_cmap3(1), label = 'Distributed circuit', linewidth = 3, zorder = 8)
plt.plot(omegas/(2*np.pi*1e9), np.abs(test_Z_equiv_LE_circuit_symbolic), color = my_cmap2(7), label = 'Equivalent circuit', linestyle = (0,(5.5,4.5)), linewidth = 3, alpha = 0.7, zorder = 11)

# plt.yscale('log')
# plt.show()
# sys.exit()

### add direct cap result

print('L_readout:', L_readout)
print('L_purcell:', L_purcell)

l_Rn = l_Rn + l_c/2
l_Rf = l_Rf + l_c/2
l_Gn = l_Gn + l_c/2
l_Gf = l_Gf + l_c/2

print('l_Rf + l_Rn:', l_Rf + l_Rn)
print('l_Gf + l_Gn:', l_Gf + l_Gn)

C_direct_val = get_eff_Cf_from_J(0, l_Gf, l_Gn, l_Rf, l_Rn, J_val, phase_vel=phase_vel, Z0=Z0)

C_direct_val = 1.38 * C_direct_val

J_cap_val = J_coupling_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, C_direct_val, phase_vel=phase_vel, Z0=Z0)

print('C_direct_val:', C_direct_val*1e15)
print('J_cap_val (MHz):', J_cap_val/(2*np.pi*1e6))

test_Z_transfer_exact = Z_transfer_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, C_direct_val, omegas, phase_vel=phase_vel, Z0=Z0)

test_Z_equiv_LE_circuit_symbolic = Z_transfer_direct_cap_equivalent_LE_circuit(l_Gf, l_Gn, l_Rf, l_Rn, C_direct_val, omegas, phase_vel=phase_vel, Z0=Z0)

# my_cmap2(7)

plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_exact), color = my_cmap3(1), label = 'cap', linewidth = 3, zorder = 5)
plt.plot(omegas/(2*np.pi*1e9), np.abs(test_Z_equiv_LE_circuit_symbolic), color = 'lightskyblue', label = 'Equiv. cap', linestyle = (0,(5.5,4.5)),  linewidth = 3, alpha = 0.7, zorder = 10)

plt.yscale('log')
# plt.legend(loc = 'upper left', fontsize = 16, frameon = False)
plt.xlabel('Frequency (GHz)', size = 20)
plt.ylabel(r'$Z_{21}$ ($\Omega$)', size = 20)
#plt.title('Z transfer function for different models')
# Customize the border settings
ax = plt.gca()
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

ax.set_xticks([8, 9, 10, 11])

ax.set_yticks([1e-2, 1, 100, 1e4])
#ax2.set_yticklabels(['0', '5', '10', '15'])

ax.tick_params(width = 2, length = 6, direction="in") ## , labelbottom=False, labelleft=False
ax.tick_params(which='minor', width = 1.5, length=3, direction='in')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

plt.xlim(omegas[0]/(2*np.pi * 1e9), omegas[-1]/(2*np.pi * 1e9))

# Adjust tick thickness and label size
# ax.tick_params(axis='both', which='both', width=2.5)
ax.tick_params(axis='both', labelsize=18)
plt.grid(True, zorder = 1)
plt.tight_layout()
plt.savefig('Z_21_circuit_models.pdf', format = 'pdf')
plt.show()
sys.exit()
