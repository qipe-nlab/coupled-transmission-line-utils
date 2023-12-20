import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
import cap_util as cap

def f(x):

    numerator = (1/x - x)**3
    denominator = (3 - 2*x**2) * (np.cos(np.pi * x/2))**2

    val = numerator/denominator

    return val

x_vals = np.linspace(0.5, 0.999, 50)

f_vals = f(x_vals)

plt.plot(x_vals, f_vals)
plt.show()

d_vals = np.linspace(1, 40, 20) * 1e-6

# print('d_vals:', d_vals)

test_Cm_vals = cap.get_Cm(d_vals)
test_Lm_vals = cap.get_Lm(d_vals)
test_Zm_vals = cap.get_Zm(d_vals)

# print('test_Cm_vals:', test_Cm_vals)
# print('test_Lm_vals:', test_Lm_vals)
# print('test_Zm_vals:', test_Zm_vals)

# # plt.plot(d_vals*1e6, test_Cm_vals * 1e15)
# # plt.show()
# # plt.plot(d_vals*1e6, test_Lm_vals * 1e9)
# # plt.show()
# # plt.plot(d_vals*1e6, test_Zm_vals)
# # plt.show()

# def RHS_eq(Zm, Zc):

#     val = Zm**2 + Zc**2/ (Zc**2 - Zm**2)

#     return val

default_phase_vel = 119919602
default_Z0 = 65

Cl = 1/(default_phase_vel * default_Z0)
Ll = default_Z0/default_phase_vel

Zc_vals = np.sqrt(Ll/(Cl + test_Cm_vals))

# RHS_test_vals = np.array([RHS_eq(Zm, Zc) for Zm, Zc in zip(test_Zm_vals, Zc_vals)])

# plt.plot(RHS_test_vals)
# plt.show()

# ###

from exact_coupled_transmission_line_eqn_solver import *

# omegas = np.linspace(2, 12, 100) * 2*np.pi * 1e9

#phase_vel=3*10**8/2.5
#Z0 = 65

phase_vel = 119919602
Z0 = 65 #65

Cl = 1/(default_phase_vel * default_Z0)
Ll = default_Z0/default_phase_vel

Cm_per_len = cap.get_Cm(15e-6)
Lm_per_len = cap.get_Lm(15e-6)

l_Rf = 1.65e-3
l_Rn = 0.75e-3
l_Gf = 2.2e-3
l_Gn = 0.21e-3
l_c = 0.45e-3 * 1

# phase_vel_c = 1/(np.sqrt(Ll*(Cl + Cm_per_len)))

# notch_val = find_notch_filter_frequency_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, Z0 = Z0, search_span = 3 * 2*np.pi*1e9, search_spacing=(10*2*np.pi*10**6))

# tau_Gf_dash = l_Gf/phase_vel + l_c/(2*phase_vel_c)
# tau_Rf_dash = l_Rf/phase_vel + l_c/(2*phase_vel_c)

# test_simple_notch_val = np.pi/(2 * (tau_Gf_dash + tau_Rf_dash))

# print('notch_val:', notch_val)
# print('test_simple_notch_val:', test_simple_notch_val)

d_vals = np.linspace(2, 25, 10) * 1e-6

print('d_vals', d_vals)

# J_vals_test = np.array([J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Lm(d), cap.get_Cm(d), phase_vel=phase_vel, Z0=Z0) for d in d_vals])

# plt.plot(d_vals, J_vals_test/(2*np.pi * 1e6))
# plt.show()

l_c_vals = np.linspace(100, 400, 15) * 1e-6

# notch_vals = np.array([find_notch_filter_frequency_analytic(l_c_val, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, Z0 = Z0, search_span = 3 * 2*np.pi*1e9, search_spacing=(10*2*np.pi*10**6))  for l_c_val in l_c_vals]) 

# plt.plot(l_c_vals, notch_vals/(2*np.pi * 1e9))
# plt.show()

cpw_length_G = l_Gf + l_Gn + l_c
cpw_length_R = l_Rf + l_Rn + l_c

omega_r = lambda_quarter_omega(cpw_length_G, phase_vel=phase_vel)
omega_p = lambda_quarter_omega(cpw_length_R, phase_vel=phase_vel)

#omega_f = find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel=phase_vel, Z0=Z0, search_span=2*2*np.pi*10**9, search_spacing=(2.5*2*np.pi*10**6))

omega_f_analytic = find_notch_filter_frequency_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel=phase_vel, Z0=Z0, search_span=2*2*np.pi*10**9, search_spacing=(2.5*2*np.pi*10**6))

omega_f_rule_of_thumb_scaled = notch_filter_frequency_rule_of_thumb(l_c, l_Gf,l_Rf, Cm_per_len, phase_vel=phase_vel, Z0=Z0, scale_phase_c = True)
omega_f_rule_of_thumb_basic = notch_filter_frequency_rule_of_thumb(l_c, l_Gf,l_Rf, Cm_per_len, phase_vel=phase_vel, Z0=Z0, scale_phase_c = False)

print('omega_r:', omega_r/(2*np.pi*1e9))
print('omega_p:', omega_p/(2*np.pi*1e9))
print('omega_f_rule_of_thumb_scaled:', omega_f_rule_of_thumb_scaled/(2*np.pi*1e9))
print('omega_f_rule_of_thumb_basic:', omega_f_rule_of_thumb_basic/(2*np.pi*1e9))

J_vals_analytic = np.array([J_coupling_analytic(l_c, l_Gf, l_Gn, l_Rf, cap.get_Cm(d_val), phase_vel=phase_vel, Z0=Z0) for d_val in d_vals])
print('J_vals_analytic:', J_vals_analytic/(2*np.pi*1e9))

####
# delta_omega = 10 * 2*np.pi
# Z_plus = Z_transfer_weak_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega_f_analytic + delta_omega/2, phase_vel=phase_vel, Z0=Z0)
# Z_minus = Z_transfer_weak_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega_f_analytic - delta_omega/2, phase_vel=phase_vel, Z0=Z0)
# Z_diff3 = (Z_plus - Z_minus)/delta_omega

# print('Z_diff:', Z_diff)
# print('Z_diff7:', Z_diff7)

# J_vals_exact_circuit = np.array([J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Lm(d_val), cap.get_Cm(d_val), phase_vel=phase_vel, Z0=Z0) for d_val in d_vals])
# J_vals_analytic = np.array([J_coupling_analytic(l_c, l_Gf, l_Gn, l_Rf, cap.get_Cm(d_val), phase_vel=phase_vel, Z0=Z0) for d_val in d_vals])

# plt.plot(d_vals * 1e6, J_vals_exact_circuit/(2*np.pi * 1e6), color = 'r')
# plt.plot(d_vals * 1e6, J_vals_analytic/(2*np.pi * 1e6), color = 'g')
# plt.show()

# plt.plot(d_vals * 1e6, 100*(J_vals_analytic-J_vals_exact_circuit)/(J_vals_exact_circuit), color = 'r')
# plt.show()

J_vals_exact_circuit = np.array([J_coupling(l_c_val, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Lm(5e-6), cap.get_Cm(5e-6), phase_vel=phase_vel, Z0=Z0) for l_c_val in l_c_vals])
J_vals_analytic = np.array([J_coupling_analytic(l_c_val, l_Gf, l_Gn, l_Rf, cap.get_Cm(5e-6), phase_vel=phase_vel, Z0=Z0) for l_c_val in l_c_vals])

plt.plot(l_c_vals, J_vals_exact_circuit/(2*np.pi * 1e6), color = 'r')
plt.plot(l_c_vals, J_vals_analytic/(2*np.pi * 1e6), color = 'g')
plt.show()

plt.plot(l_c_vals, 100*(J_vals_analytic-J_vals_exact_circuit)/(J_vals_exact_circuit), color = 'r')
plt.show()
###