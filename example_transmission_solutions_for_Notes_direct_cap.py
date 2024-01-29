import seaborn as sns
from scipy import constants

hex_codes = sns.color_palette().as_hex()
hex_codes2 = sns.color_palette('husl', 9).as_hex()

from matplotlib.colors import ListedColormap

#flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
my_cmap = ListedColormap(hex_codes)
my_cmap2 = ListedColormap(hex_codes2)
my_cmap3 = ListedColormap(["#9b59b6", "#34495e"])

### example solutions

from exact_coupled_transmission_line_eqn_solver import *
import cap_util as cap

### param set 1: for high freq range

phase_vel = 1.2e8 #constants.c / np.sqrt(6.351)
Z0 = 65

d_val = 10e-6 # 12.5*1e-6 # separation between the coupled sections

CJ = 2e-15

l_Rf = 1.7e-3
l_Rn = 1.21e-3
l_Gf = 1.2e-3
l_Gn = 1.75e-3

test_J_value = J_coupling_direct_cap_analytic(l_Gf, l_Gn, l_Rf, l_Rn, CJ, phase_vel=phase_vel, Z0=Z0)

print('test_J_value:',test_J_value/(2*np.pi*1e6))

L_readout = l_Gn + l_Gf

readout_omega = lambda_quarter_omega(L_readout, phase_vel=phase_vel)

L_purcell = l_Rn  + l_Rf

purcell_omega = lambda_quarter_omega(L_purcell, phase_vel=phase_vel)

print('readout_omega:', readout_omega / (2*np.pi*1e9))
print('purcell_omega:', purcell_omega / (2*np.pi*1e9))

print('l_Rf:', l_Rf)
print('l_Rn:', l_Rn)
print('l_Gf:', l_Gf)
print('l_Gn:', l_Gn)

test_omega_1 = lambda_quarter_omega(l_Gf + l_Gn, phase_vel=phase_vel)

test_omega_2 = lambda_quarter_omega(l_Rf + l_Rn, phase_vel=phase_vel)

# test_J_val = J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm)

# print('test_notch_freq (GHz):', test_notch_freq/(2*np.pi*1e9))
print('test_omega_1 (GHz):', test_omega_1/(2*np.pi*1e9))
print('test_omega_2 (GHz):', test_omega_2/(2*np.pi*1e9))
# print('test_J_val (MHz):', test_J_val/(2*np.pi*1e6))

omegas = np.arange(7.5, 11, 0.02) * 1e9 * 2*np.pi

Z11_exact = Z_input_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, CJ, omegas, phase_vel=phase_vel, Z0=Z0)
Z11_equivalent_circuit_symbolic = Z_input_direct_cap_equivalent_LE_circuit(l_Gf, l_Gn, l_Rf, l_Rn, CJ, omegas, phase_vel=phase_vel, Z0=Z0)

plt.plot(omegas/(2*np.pi * 1e9), np.imag(Z11_exact), color = my_cmap3(1), label = 'Exact solution', linewidth = 3)
plt.plot(omegas/(2*np.pi * 1e9), np.imag(Z11_equivalent_circuit_symbolic), color = my_cmap2(7), label = 'Equivalent circuit - symbolic', linewidth = 3, alpha = 0.6)

plt.yscale('log')

plt.xlabel('Frequency (GHz)', size = 20)
plt.ylabel(r'$Z_{11}$ ($\Omega$)' + r' $(\Omega)$', size = 20)
plt.legend(fontsize = 16)
#plt.title('Radiative T1 limit through detector line')

ax = plt.gca()

ax.tick_params(axis='both', which='major', labelsize=16)

# ax.set_yticks([1e-6,1e-3,1,1e3]) 
# ax.set_yticklabels(['1 ns','1 us','1 ms', '1 s'])
# ax.set_ylim([0.5e-6, 5e3])

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)

ax.tick_params(width=3)
plt.tight_layout()
plt.show()

###

test_Z_transfer_exact = Z_transfer_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, CJ, omegas, phase_vel=phase_vel, Z0=Z0)

test_Z_equiv_LE_circuit_symbolic = Z_transfer_direct_cap_equivalent_LE_circuit(l_Gf, l_Gn, l_Rf, l_Rn, CJ, omegas, phase_vel=phase_vel, Z0=Z0)

# my_cmap2(7)

plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_exact), color = my_cmap3(1), label = 'Distributed circuit', linewidth = 3)
#plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_weak_coupling), color = 'r', linestyle = '--', label = 'Z transfer weak coupling model')
#plt.plot(omegas/(2*np.pi*1e9), np.abs(test_Z_equiv_LE_circuit), color = my_cmap2(7), label = 'Equivalent circuit', linewidth = 3, alpha = 0.85)
plt.plot(omegas/(2*np.pi*1e9), np.abs(test_Z_equiv_LE_circuit_symbolic), color = my_cmap2(7), label = 'Equivalent circuit - sybc.', linewidth = 3, alpha = 0.5)

#plt.vlines(test_notch_freq_analytic/(2*np.pi*1e9), 0.008, 2, color = my_cmap(7), linestyle = 'dotted', linewidth = 3)
#plt.vlines(test_notch_freq_rule_of_thumb/(2*np.pi*1e9), 0.008, 2, color = my_cmap(7), linestyle = 'dotted', linewidth = 3, label = 'Notch eq.',zorder=10)

plt.yscale('log')
plt.legend(loc = 'upper left', fontsize = 14)
plt.xlabel('Frequency (GHz)', size = 20)
plt.ylabel(r'$Z_{21}$ ($\Omega$)', size = 20)
#plt.title('Z transfer function for different models')
# Customize the border settings
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(3)
ax.spines['left'].set_linewidth(3)

ax.set_xticks([8, 9, 10, 11])

ax.set_yticks([1e-2, 1, 100, 1e4])
#ax2.set_yticklabels(['0', '5', '10', '15'])

# Adjust tick thickness and label size
ax.tick_params(axis='both', which='both', width=2.5)
ax.tick_params(axis='both', labelsize=20)
plt.tight_layout()

plt.show()

#sys.exit()

C_q = 50e-15

### param set 1: for high freq range
C_g = 5.5e-15 #3.65e-15
C_ext = 10e-15 #24e-15

### param set 2: for low freq range
# C_g = 6e-15
# C_ext = 52.5e-15

omega_q = test_notch_freq 

test_g_val = g_coupling(omega_q, C_q, C_g, l_c, l_Gf, l_Gn, phase_vel=3*10**8/2.5, Z0=65)
test_J_val = J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=3*10**8/2.5, Z0=65)
test_J_val_symbolic = J_coupling_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65)

test_k_val = k_readout(C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)

print('test_g_val (MHz):', test_g_val/(2*np.pi*1e6))
print('test_J_val (MHz):', test_J_val/(2*np.pi*1e6))
print('test_J_val_symbolic (MHz):', test_J_val_symbolic/(2*np.pi*1e6))

J_ratios = test_J_val/test_J_val_symbolic

print('J_ratios:', J_ratios)
print('Z_ratios:', Z_ratios)

print('test_k_val (MHz):', test_k_val/(2*np.pi*1e6))

omegas = np.arange(7.5, 11, 0.02) * 1e9 * 2*np.pi

test_T1_radiative_exact = qubit_radiative_decay_sym_3_lines_exact(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)
test_T1_radiative_equivalent_LE_circuit = qubit_radiative_decay_equivalent_LE_circuit(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)
test_T1_radiative_equivalent_LE_circuit_without_notch = qubit_radiative_decay_equivalent_LE_circuit_without_notch(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)
test_T1_radiative_equivalent_LE_circuit_single_resonator = qubit_radiative_decay_equivalent_LE_circuit_single_resonator(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)

plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_exact * 1e3, color = my_cmap3(1), linewidth = 3, label = 'with intrinsic notch - exact')
plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit * 1e3, color = my_cmap2(7), linewidth = 3, label = 'with intrinsic notch - e.c.', alpha = 0.5)
plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit_without_notch * 1e3, color = my_cmap2(7), linewidth = 3, linestyle = '--', label = 'without intrinsic notch - e.c.', alpha = 0.5)
# plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit_single_resonator * 1e3, color = 'b', linestyle = '--', label = 'equiv. single resonator circuit')

plt.yscale('log')
plt.xlabel('Frequency (GHz)', size = 20)
plt.ylabel(r'Purcell decay $T_1$ limit', size = 20)
plt.legend(fontsize = 16)
#plt.title('Radiative T1 limit through detector line')

ax = plt.gca()

ax.tick_params(axis='both', which='major', labelsize=16)

ax.set_yticks([1e-6,1e-3,1,1e3]) 
ax.set_yticklabels(['1 ns','1 us','1 ms', '1 s'])
ax.set_ylim([0.5e-6, 5e3])

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)

ax.tick_params(width=3)
plt.tight_layout()

plt.show()

enhancement_factor = test_T1_radiative_exact/test_T1_radiative_equivalent_LE_circuit_without_notch
enhancement_factor_approx = test_T1_radiative_equivalent_LE_circuit/test_T1_radiative_equivalent_LE_circuit_without_notch

test_enhancement_bandwidth = notch_enhancement_bandwidth(l_c, l_Gf, l_Gn, l_Rf, l_Rn, T1_enhancement_fact = 10, phase_vel = 3*10**8/2.5)

print('test_enhancement_bandwidth (GHz):', test_enhancement_bandwidth / (2*np.pi*1e9))

plt.plot(omegas/(2*np.pi * 1e9), enhancement_factor, color = my_cmap3(1), linewidth = 3, label = 'exact soltution')
plt.plot(omegas/(2*np.pi * 1e9), enhancement_factor_approx, color = my_cmap2(7), linewidth = 3, label = 'equivalent circuit solution', alpha = 0.5)

omega_r = lambda_quarter_omega(l_c + l_Gf + l_Gn, phase_vel=phase_vel)
omega_p = lambda_quarter_omega(l_c + l_Rf + l_Rn, phase_vel=phase_vel)

predicted_notch_enhancement = enhancement_factor_symbolic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, omegas)
#plt.plot(omegas/(2*np.pi * 1e9), predicted_notch_enhancement, color = 'g')

exact_bandwidth = find_notch_enhancement_bandwidth(omegas, enhancement_factor, 10)

print('exact_bandwidth (GHz):', exact_bandwidth / (2*np.pi*1e9))

plt.yscale('log')
plt.xlabel('Frequency (GHz)', size = 20)
plt.ylabel(r'Purcell $T_1$ enhancement', size = 20)
plt.legend(fontsize = 16)

ax = plt.gca()

ax.tick_params(axis='both', which='major', labelsize=16)

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)

ax.tick_params(width=3)
plt.tight_layout()

plt.hlines(10, 6, 11, color = 'k', linestyle = '--')
plt.yscale('log')
plt.show()