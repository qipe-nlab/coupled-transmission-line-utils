### example solutions

from exact_coupled_transmission_line_eqn_solver import *
import cap_util as cap

### param set 1: for high freq range

phase_vel = 3*10**8/2.5
Z0 = 65

d_val = 10e-6 # 12.5*1e-6 # separation between the coupled sections

Cm = cap.get_Cm(d_val)
Lm = cap.get_Lm(d_val)

l_Rf = 1.5e-3
l_Rn = 1.175e-3
l_Gf = 1.925e-3
l_Gn = 0.75e-3
l_c = 0.3e-3

L_readout = l_Gn + l_c + l_Gf

readout_omega = lambda_quarter_omega(L_readout, phase_vel=phase_vel)

L_purcell = l_Rn + l_c + l_Rf

purcell_omega = lambda_quarter_omega(L_purcell, phase_vel=phase_vel)

print('readout_omega:', readout_omega / (2*np.pi*1e9))
print('purcell_omega:', purcell_omega / (2*np.pi*1e9))

# omega = 8 * 2*np.pi*1e9
#C = 100 * 1e-15

# L = 1/(C*omega**2)

# print('L:', L)

# L = 10 * 1e-9

# C = 1/(L*omega**2)

# print('C:', C)

### param set 2: for low freq range

# scale_fac = 1.6

# Lm = 1.31e-9*1e1*2.2
# Cm = 5e-15*1e3*2.2
# l_Rf = 1.15e-3 * scale_fac
# l_Rn = 1.45e-3 * scale_fac
# l_Gf = 1.7e-3 * scale_fac
# l_Gn = 0.9e-3 * scale_fac
# l_c = 0.5e-3 * scale_fac

print('Lm:', Lm)
print('Cm:', Cm)
print('l_Rf:', l_Rf)
print('l_Rn:', l_Rn)
print('l_Gf:', l_Gf)
print('l_Gn:', l_Gn)
print('l_c:', l_c)
#sys.exit()

test_notch_freq = find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=phase_vel, Z0=Z0)

test_notch_freq_analytic = find_notch_filter_frequency_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=phase_vel, Z0=Z0)

test_notch_freq_rule_of_thumb = notch_filter_frequency_rule_of_thumb(l_c, l_Gf,l_Rf, Cm, phase_vel=phase_vel, Z0=Z0)

# notch_freq_agreement_percent = 100*(test_notch_freq_analytic - test_notch_freq )/ test_notch_freq

# print('notch_freq_agreement_percent:', notch_freq_agreement_percent)

# testing LE coupling between TLs
# Cc = 250*1e-15
# Lc = 3e-9

# omegas = np.arange(2, 10, 0.02) * 1e9 * 2*np.pi

# l_Gn = 2.5e-3 * scale_fac
# l_Gf = 0.4e-3 * scale_fac

# l_Rn = 2.5e-3 * scale_fac
# l_Rf = 0.4e-3 * scale_fac

# C_gval, L_gval = lumped_model_Cg_and_Lg_LE_coupling(3*10**8/2.5, 65, l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc)

# print('C_gval, L_gval:', C_gval, L_gval)

# lumped_element_values = get_lumped_elements_LE_coupling(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, phase_vel=3*10**8/2.5, Z0=65)

# test_res_coupling = lumped_model_resonator_coupling(*lumped_element_values)

# print('lumped_element_values:', lumped_element_values)

# print('test_res_coupling (MHz):', test_res_coupling/(2*np.pi*1e6))

# test_Zvals = Z_transfer_LE_coupled_lines(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, omegas, phase_vel=3*10**8/2.5, Z0=65)

# test_Zequiv_vals = Z_transfer_equivalent_LE_circuit(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, omegas, phase_vel=3*10**8/2.5, Z0=65)

# plt.plot(omegas/(2*np.pi*1e9), np.abs(np.imag(test_Zvals)))
# plt.plot(omegas/(2*np.pi*1e9), np.abs(np.imag(test_Zequiv_vals)), linestyle = '--')

# plt.yscale('log')
# plt.show()
# print('test_Zvals:', test_Zvals)

# #sys.exit()

test_omega_1 = lambda_quarter_omega(l_c + l_Gf + l_Gn, phase_vel=phase_vel)

test_omega_2 = lambda_quarter_omega(l_c + l_Rf + l_Rn, phase_vel=phase_vel)

# test_J_val = J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm)

# print('test_notch_freq (GHz):', test_notch_freq/(2*np.pi*1e9))
print('test_omega_1 (GHz):', test_omega_1/(2*np.pi*1e9))
print('test_omega_2 (GHz):', test_omega_2/(2*np.pi*1e9))
# print('test_J_val (MHz):', test_J_val/(2*np.pi*1e6))

omegas = np.arange(6.5, 11, 0.02) * 1e9 * 2*np.pi

test_Z_transfer_exact = Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

test_S_transfer_sym_3_lines_exact = S_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

print('test_S_transfer_sym_3_lines_exact:', test_S_transfer_sym_3_lines_exact)

S_transfer = test_S_transfer_sym_3_lines_exact[:, 0, 1]

print('S_transfer:', S_transfer)

plt.plot(omegas/(2*np.pi*1e9), np.abs(S_transfer))
plt.yscale('log')
plt.show()

test_Z_transfer_weak_coupling = Z_transfer_weak_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

test_Z_equiv_LE_circuit = Z_transfer_equivalent_LE_circuit(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

# print('test_Z_transfer_exact:', test_Z_transfer_exact[:10])
# print('test_Z_transfer_weak_coupling:', test_Z_transfer_weak_coupling[:10])
# print('test_Z_transfer_equiv_LE_circuit:', test_Z_equiv_LE_circuit[:10])

import seaborn as sns
hex_codes = sns.color_palette().as_hex()
hex_codes2 = sns.color_palette('husl', 9).as_hex()

from matplotlib.colors import ListedColormap

#flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
my_cmap = ListedColormap(hex_codes)
my_cmap2 = ListedColormap(hex_codes2)
my_cmap3 = ListedColormap(["#9b59b6", "#34495e"])

# my_cmap2(7)

plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_exact), color = my_cmap3(1), label = 'Distributed circuit', linewidth = 3)
#plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_weak_coupling), color = 'r', linestyle = '--', label = 'Z transfer weak coupling model')
plt.plot(omegas/(2*np.pi*1e9), np.abs(test_Z_equiv_LE_circuit), color = my_cmap2(7), label = 'Equivalent circuit', linewidth = 3, alpha = 0.85)

#plt.vlines(test_notch_freq_analytic/(2*np.pi*1e9), 0.008, 2, color = my_cmap(7), linestyle = 'dotted', linewidth = 3)
plt.vlines(test_notch_freq_rule_of_thumb/(2*np.pi*1e9), 0.008, 2, color = my_cmap(7), linestyle = 'dotted', linewidth = 3, label = 'Notch eq. (x)',zorder=10)

plt.yscale('log')
plt.legend(loc = 'lower right', fontsize = 16)
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
test_k_val = k_readout(C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)

print('test_g_val (MHz):', test_g_val/(2*np.pi*1e6))
print('test_J_val (MHz):', test_J_val/(2*np.pi*1e6))
print('test_k_val (MHz):', test_k_val/(2*np.pi*1e6))

omegas = np.arange(6, 11, 0.02) * 1e9 * 2*np.pi

test_T1_radiative_exact = qubit_radiative_decay_sym_3_lines_exact(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)
test_T1_radiative_equivalent_LE_circuit = qubit_radiative_decay_equivalent_LE_circuit(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)
test_T1_radiative_equivalent_LE_circuit_without_notch = qubit_radiative_decay_equivalent_LE_circuit_without_notch(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)
test_T1_radiative_equivalent_LE_circuit_single_resonator = qubit_radiative_decay_equivalent_LE_circuit_single_resonator(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)

plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_exact * 1e3, color = 'k', label = 'with intrinsic notch')
plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit * 1e3, color = 'r', linestyle = '--')
plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit_without_notch * 1e3, color = 'r', linestyle = '--', label = 'without intrinsic notch')
# plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit_single_resonator * 1e3, color = 'b', linestyle = '--', label = 'equiv. single resonator circuit')

plt.yscale('log')
plt.xlabel('Frequency (GHz)', size = 20)
plt.ylabel(r'Purcell decay $T_1$ limit', size = 20)
plt.legend()
#plt.title('Radiative T1 limit through detector line')

ax = plt.gca()

ax.tick_params(axis='both', which='major', labelsize=15)

ax.set_yticks([1e-6,1e-3,1,1e3]) 
ax.set_yticklabels(['1 ns','1 us','1 ms', '1 s'])
ax.set_ylim([0.5e-6, 5e3])

# plt.yticks([1e-5,1e-4,1e-2,1e-1,1e1, 1e2],
#            ['',"", "", "", '', ''], minor=True)

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)

ax.tick_params(width=3)
plt.tight_layout()

plt.show()

enhancement_factor = test_T1_radiative_exact/test_T1_radiative_equivalent_LE_circuit_without_notch
enhancement_factor_approx = test_T1_radiative_equivalent_LE_circuit/test_T1_radiative_equivalent_LE_circuit_without_notch

test_enhancement_factor = notch_enhancement_bandwidth(l_c, l_Gf, l_Gn, l_Rf, l_Rn, T1_enhancement_fact = 10, phase_vel = 3*10**8/2.5)

print('test_enhancement_factor (GHz):', test_enhancement_factor / (2*np.pi*1e9))

plt.plot(omegas/(2*np.pi * 1e9), enhancement_factor, color = 'r')
plt.plot(omegas/(2*np.pi * 1e9), enhancement_factor_approx, color = 'b')

omega_r = lambda_quarter_omega(l_c + l_Gf + l_Gn, phase_vel=phase_vel)
omega_p = lambda_quarter_omega(l_c + l_Rf + l_Rn, phase_vel=phase_vel)

predicted_notch_enhancement = enhancement_factor_symbolic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, omegas)
plt.plot(omegas/(2*np.pi * 1e9), predicted_notch_enhancement, color = 'g')


plt.hlines(10, 6, 11, color = 'k', linestyle = '--')
plt.yscale('log')
plt.show()


omega_n = 8 
omega_r = 10
omega_p = 10

test = notch_enhancement_bandwidth_by_omegas(omega_n, omega_r, omega_p, T1_enhancement_fact = 10)

print(test)