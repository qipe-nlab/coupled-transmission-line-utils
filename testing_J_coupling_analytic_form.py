import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
import cap_util as cap

import seaborn as sns
from scipy import constants

hex_codes = sns.color_palette().as_hex()
hex_codes2 = sns.color_palette('husl', 9).as_hex()

from matplotlib.colors import ListedColormap

#flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
my_cmap = ListedColormap(hex_codes)
my_cmap2 = ListedColormap(hex_codes2)
my_cmap3 = ListedColormap(["#9b59b6", "#34495e"])

d_vals = np.linspace(1, 40, 20) * 1e-6

# print('d_vals:', d_vals)

test_Cm_vals = cap.get_Cm(d_vals)
test_Lm_vals = cap.get_Lm(d_vals)
test_Zm_vals = cap.get_Zm(d_vals)

# default_phase_vel = 119919602
# default_Z0 = 65

from exact_coupled_transmission_line_eqn_solver import *

phase_vel = 1.2e8
Z0 = 65

Cm_per_len = cap.get_Cm(10e-6)
Lm_per_len = cap.get_Lm(10e-6)

l_Rf = 1.7e-3
l_Rn = 0.96e-3
l_Gf = 1.4e-3
l_Gn = 1.3e-3
l_c = 0.25e-3

# l_Rf = 1.85e-3
# l_Rn = 0.75e-3 + 0.35e-3
# l_Gf = 2.4e-3 
# l_Gn = 0.21e-3 + 0.15e-3
# l_c = 0.35e-3 

# phase_vel_c = 1/(np.sqrt(Ll*(Cl + Cm_per_len)))

# notch_val = find_notch_filter_frequency_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, Z0 = Z0, search_span = 3 * 2*np.pi*1e9, search_spacing=(10*2*np.pi*10**6))

# tau_Gf_dash = l_Gf/phase_vel + l_c/(2*phase_vel_c)
# tau_Rf_dash = l_Rf/phase_vel + l_c/(2*phase_vel_c)

# test_simple_notch_val = np.pi/(2 * (tau_Gf_dash + tau_Rf_dash))

# print('notch_val:', notch_val)
# print('test_simple_notch_val:', test_simple_notch_val)

d_vals = np.linspace(2, 10, 20) * 1e-6

print('d_vals', d_vals)

l_c_vals = np.linspace(100, 400, 15) * 1e-6

cpw_length_G = l_Gf + l_Gn + l_c
cpw_length_R = l_Rf + l_Rn + l_c

omega_r = lambda_quarter_omega(cpw_length_G, phase_vel=phase_vel)
omega_p = lambda_quarter_omega(cpw_length_R, phase_vel=phase_vel)

#omega_f_analytic = find_notch_filter_frequency_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel=phase_vel, Z0=Z0, search_span=2*2*np.pi*10**9, search_spacing=(2.5*2*np.pi*10**6))

omega_f_rule_of_thumb_scaled = notch_filter_frequency_rule_of_thumb(l_c, l_Gf,l_Rf, Cm_per_len, phase_vel=phase_vel, Z0=Z0)
omega_f_rule_of_thumb_basic = notch_filter_frequency_rule_of_thumb(l_c, l_Gf, l_Rf, phase_vel=phase_vel, Z0=Z0)

print('omega_r:', omega_r/(2*np.pi*1e9))
print('omega_p:', omega_p/(2*np.pi*1e9))
print('omega_f_rule_of_thumb_scaled:', omega_f_rule_of_thumb_scaled/(2*np.pi*1e9))
print('omega_f_rule_of_thumb_basic:', omega_f_rule_of_thumb_basic/(2*np.pi*1e9))

J_vals_solve_eq_c = np.array([J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Lm(d_val), cap.get_Cm(d_val), phase_vel=phase_vel, Z0=Z0) for d_val in d_vals])
print('J_vals_solve_eq_c:', J_vals_solve_eq_c/(2*np.pi*1e9))

J_vals_analytic_raw = np.array([J_coupling_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Cm(d_val), phase_vel=phase_vel, Z0=Z0, simplified=False) for d_val in d_vals])
print('J_vals_analytic_raw:', J_vals_analytic_raw/(2*np.pi*1e9))

J_vals_analytic_simplified = np.array([J_coupling_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Cm(d_val), phase_vel=phase_vel, Z0=Z0, simplified=True) for d_val in d_vals])
print('J_vals_analytic_simplified:', J_vals_analytic_simplified/(2*np.pi*1e9))

ax = plt.gca()

ax.tick_params(axis='both', which='major', labelsize=16)

# ax.set_yticks([1e-6,1e-3,1,1e3]) 
# ax.set_yticklabels(['1 ns','1 us','1 ms', '1 s'])
# ax.set_ylim([0.5e-6, 5e3])

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)

ax.tick_params(width=3)

plt.xlabel(r'$d$ (um)', size = 20)
plt.ylabel(r'$J/2\pi$ (MHz)', size = 20)

plt.plot(d_vals * 1e6, J_vals_solve_eq_c/(2*np.pi*1e6), linestyle = 'none', marker = 'o', color = my_cmap3(1), label = 'Numeric', linewidth = 3, markersize = 9)
#plt.plot(d_vals * 1e6, J_vals_analytic_raw/(2*np.pi*1e6), color = 'b')
plt.plot(d_vals * 1e6, J_vals_analytic_simplified/(2*np.pi*1e6), linestyle = 'none', marker = 'o', color = my_cmap2(7), label = 'Symbolic', linewidth = 3, markersize = 9)

plt.legend(loc = 'upper right', fontsize = 14)
plt.tight_layout()
plt.show()

sys.exit()

l_c_exp = 340e-6 
s_exp = 5e-6
## SW Q0:
# omega_n_exp = 8.270*2*np.pi*1e9
# omega_r_exp = 10.264*2*np.pi*1e9
## SW Q1:
# omega_n_exp = 8.960*2*np.pi*1e9
# omega_r_exp = 10.672*2*np.pi*1e9
## SW Q2:
omega_n_exp = 8.690*2*np.pi*1e9
omega_r_exp = 10.479*2*np.pi*1e9
## SW Q3:
# omega_n_exp = 8.270*2*np.pi*1e9
# omega_r_exp = 10.04*2*np.pi*1e9

## CD36 device design

## U1 Q1
s_exp = 3.8e-6
l_c_exp = 330e-6 
omega_n_exp = 9.1*2*np.pi*1e9
omega_r_exp = 10.650*2*np.pi*1e9
omega_p_exp = 10.650*2*np.pi*1e9

Cm_per_len_exp = cap.get_Cm(s_exp)

exp_J_analytic_prediction = J_coupling_analytic_by_freqs(omega_r_exp, omega_p_exp, omega_n_exp, l_c_exp, Cm_per_len_exp, phase_vel=phase_vel, Z0=Z0)

print('exp_J_analytic_prediction (MHz):', exp_J_analytic_prediction / (2*np.pi*1e6))

#sys.exit()

####
# delta_omega = 10 * 2*np.pi
# Z_plus = Z_transfer_weak_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega_f_analytic + delta_omega/2, phase_vel=phase_vel, Z0=Z0)
# Z_minus = Z_transfer_weak_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega_f_analytic - delta_omega/2, phase_vel=phase_vel, Z0=Z0)
# Z_diff3 = (Z_plus - Z_minus)/delta_omega

# print('Z_diff:', Z_diff)
# print('Z_diff7:', Z_diff7)

J_vals_exact_circuit = np.array([J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Lm(d_val), cap.get_Cm(d_val), phase_vel=phase_vel, Z0=Z0) for d_val in d_vals])
J_vals_analytic = np.array([J_coupling_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Cm(d_val), phase_vel=phase_vel, Z0=Z0) for d_val in d_vals])

plt.plot(d_vals * 1e6, J_vals_exact_circuit/(2*np.pi * 1e6), color = 'r')
plt.plot(d_vals * 1e6, J_vals_analytic/(2*np.pi * 1e6), color = 'g')
plt.show()

plt.plot(d_vals * 1e6, 100*(J_vals_analytic-J_vals_exact_circuit)/(J_vals_exact_circuit), color = 'r')
plt.show()

J_vals_exact_circuit = np.array([J_coupling(l_c_val, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Lm(10e-6), cap.get_Cm(10e-6), phase_vel=phase_vel, Z0=Z0) for l_c_val in l_c_vals])
J_vals_analytic = np.array([J_coupling_analytic(l_c_val, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Cm(10e-6), phase_vel=phase_vel, Z0=Z0) for l_c_val in l_c_vals])

plt.plot(l_c_vals, J_vals_exact_circuit/(2*np.pi * 1e6), color = 'r')
plt.plot(l_c_vals, J_vals_analytic/(2*np.pi * 1e6), color = 'g')
plt.show()

plt.plot(l_c_vals, 100*(J_vals_analytic-J_vals_exact_circuit)/(J_vals_exact_circuit), color = 'r')
plt.show()
###
