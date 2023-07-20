### example solutions

from exact_coupled_transmission_line_eqn_solver import *

### param set 1: for high freq range

phase_vel = 3*10**8/2.5
Z0 = 65

Lm = 1.31e-9*1e1*1.8
Cm = 5e-15*1e3*1.8
l_Rf = 1.2e-3
l_Rn = 1.2e-3
l_Gf = 1.2e-3
l_Gn = 1.2e-3
l_c = 0.5e-3

### param set 2: for low freq range

# scale_fac = 1.6

# Lm = 1.31e-9*1e1*2.2
# Cm = 5e-15*1e3*2.2
# l_Rf = 1.15e-3 * scale_fac
# l_Rn = 1.45e-3 * scale_fac
# l_Gf = 1.7e-3 * scale_fac
# l_Gn = 0.9e-3 * scale_fac
# l_c = 0.5e-3 * scale_fac

test_notch_freq = find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=3*10**8/2.5, Z0=65)

test_omega_1 = lambda_quarter_omega(l_c + l_Gf + l_Gn, phase_vel=3*10**8/2.5)

test_omega_2 = lambda_quarter_omega(l_c + l_Rf + l_Rn, phase_vel=3*10**8/2.5)

test_J_val = J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm)

print('test_notch_freq (GHz):', test_notch_freq/(2*np.pi*1e9))
print('test_omega_1 (GHz):', test_omega_1/(2*np.pi*1e9))
print('test_omega_2 (GHz):', test_omega_2/(2*np.pi*1e9))
print('test_J_val (MHz):', test_J_val/(2*np.pi*1e6))

omegas = np.arange(6, 11, 0.02) * 1e9 * 2*np.pi

test_Z_transfer_exact = Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65)

test_Z_transfer_weak_coupling = Z_transfer_weak_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65)

test_Z_equiv_LE_circuit = Z_transfer_equivalent_LE_circuit(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65)

# print('test_Z_transfer_exact:', test_Z_transfer_exact[:10])
# print('test_Z_transfer_weak_coupling:', test_Z_transfer_weak_coupling[:10])
# print('test_Z_transfer_equiv_LE_circuit:', test_Z_equiv_LE_circuit[:10])

plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_exact), color = 'k', label = 'Z transfer exact multiconductor TL model')
plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_weak_coupling), color = 'r', linestyle = '--', label = 'Z transfer weak coupling model')
plt.plot(omegas/(2*np.pi*1e9), np.abs(test_Z_equiv_LE_circuit), color = 'b', linestyle = '-.', label = 'Z transfer equivalent circuit model')
plt.yscale('log')
plt.legend()
plt.xlabel('frequency (GHz)')
plt.ylabel('Z transfer (Ohm)')
plt.title('Z transfer function for different models')
plt.show()

C_q = 50e-15

### param set 1: for high freq range
C_g = 5e-15 #3.65e-15
C_ext = 23.5e-15

### param set 2: for low freq range
# C_g = 8e-15
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

plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_exact * 1e3, color = 'k', label = 'solution with inductive-capacitive coupling')
#plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit * 1e6, color = 'r', linestyle = '--')
plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit_without_notch * 1e3, color = 'r', linestyle = '--', label = 'solution w/o inductive-capacitive coupling')
# plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit_single_resonator * 1e3, color = 'b', linestyle = '--', label = 'equiv. single resonator circuit')

plt.yscale('log')
plt.xlabel('Frequency (GHz)')
plt.ylabel(r'$T_1$ (ms)')
plt.legend()
plt.title('Radiative T1 limit through detector line')
plt.show()