import numpy as np
import scipy as sp
from scipy import constants
import matplotlib.pyplot as plt

import sys
import cap_util as cap

d_vals = np.linspace(1, 40, 20) * 1e-6

# print('d_vals:', d_vals)

Cm_vals = cap.get_Cm(d_vals)
Lm_vals = cap.get_Lm(d_vals)
Zm_vals = cap.get_Zm(d_vals)

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

default_phase_vel = constants.c / np.sqrt(6.351)
default_Z0 = 65.62

Cl = 1/(default_phase_vel * default_Z0)
Ll = default_Z0/default_phase_vel

Zc_vals = np.sqrt(Ll/(Cl + Cm_vals))

plt.plot(d_vals* 1e6, Zc_vals, color = 'r', label = r'$Z_c$', marker = 'o', markersize = 9, linestyle = 'None')
plt.plot(d_vals* 1e6, Zm_vals, color = 'purple', label = r'$Z_m$', marker = 'o', markersize = 9, linestyle = 'None')
# plt.xlabel(r'$d$ (um)')
# plt.ylabel(r'$Z$ ($\Omega$)',)
#plt.legend(loc = 'lower right')

ax = plt.gca()

# Customize the border settings
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2.5)
ax.spines['left'].set_linewidth(2.5)

# Adjust tick thickness and label size
ax.tick_params(axis='both', which='both', width=2.5)
ax.tick_params(axis='both', labelsize=25)

ax.set_xticks([0, 5, 10, 15, 20])
ax.set_xticklabels(['0', '5', '10', '15', '20'])

# Set axis labels and font size
ax.set_xlabel(r'$d$ (um)', fontsize=30)
ax.set_ylabel(r'$Z$ ($\Omega$)', fontsize=30)

# Create a second y-axis for the new dataset
ax2 = ax.twinx()
line_c, = ax2.plot(d_vals* 1e6, Cm_vals * 1e12, label=r'$C_m$', marker = '>', markersize = 9,  color='g', linestyle = 'None')
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['right'].set_linewidth(2.5)
ax2.tick_params(axis='both', which='both', width=2.5)
ax2.tick_params(axis='both', labelsize=25)
ax2.set_ylabel(r'$C_m$ (pF)', fontsize=30)
ax2.set_ylim(0, 16)

ax2.set_yticks([0, 5, 10, 15])
ax2.set_yticklabels(['0', '5', '10', '15'])


# Add legend with fontsize 20
ax.legend(loc = 'upper right', fontsize=25)
ax.set_xlim(0, 20)
ax.set_ylim(60, 68)
plt.tight_layout()
plt.show()


print('Cl:', Cl)
print('Ll:', Ll)

coupling_coefficient = cap.get_Lm(d_vals)/Ll

plt.plot(d_vals* 1e6, coupling_coefficient)
plt.show()

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

Cm_per_len = cap.get_Cm(10e-6)
Lm_per_len = cap.get_Lm(10e-6)

l_Rf = 1.15e-3 
l_Rn = 1.35e-3 
l_Gf = 1.7e-3 
l_Gn = 0.925e-3 
l_c = 0.325e-3 

# test_vals = voltage_transmission_coupled_lines_debug(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas)

phase_vel_c = 1/(np.sqrt(Ll*(Cl + Cm_per_len)))

notch_val = find_notch_filter_frequency_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, Z0 = Z0, search_span = 3 * 2*np.pi*1e9, search_spacing=(10*2*np.pi*10**6))

tau_Gf_dash = l_Gf/phase_vel + l_c/(2*phase_vel_c)
tau_Rf_dash = l_Rf/phase_vel + l_c/(2*phase_vel_c)

test_simple_notch_val = np.pi/(2 * (tau_Gf_dash + tau_Rf_dash))

print('notch_val:', notch_val)
print('test_simple_notch_val:', test_simple_notch_val)

d_vals = np.linspace(2, 25, 10) * 1e-6

print('d_vals', d_vals)
#sys.exit()

# J_vals_test = np.array([J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Lm(d), cap.get_Cm(d), phase_vel=phase_vel, Z0=Z0) for d in d_vals])

# plt.plot(d_vals, J_vals_test/(2*np.pi * 1e6))
# plt.show()

l_c_vals = np.linspace(100, 500, 10) * 1e-6

# notch_vals = np.array([find_notch_filter_frequency_analytic(l_c_val, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, Z0 = Z0, search_span = 3 * 2*np.pi*1e9, search_spacing=(10*2*np.pi*10**6))  for l_c_val in l_c_vals]) 

# plt.plot(l_c_vals, notch_vals/(2*np.pi * 1e9))
# plt.show()

J_vals_test = np.array([J_coupling(l_c*1.5, l_Gf, l_Gn, l_Rf, l_Rn, cap.get_Lm(d_val), cap.get_Cm(d_val), phase_vel=phase_vel, Z0=Z0) for d_val in d_vals])

plt.plot(l_c_vals, J_vals_test/(2*np.pi * 1e6))
plt.show()

sys.exit()

# plt.plot(omegas/(2*np.pi*1e9), np.real(test_vals))
# plt.plot(omegas/(2*np.pi*1e9), np.imag(test_vals))
# #plt.vlines(notch_val/(2*np.pi*1e9), 0, 0.01)
# plt.show()

# print('test_vals:', test_vals)

###

## quick test

def test_func(a,b):

    val = (1 + np.exp(-1j*(a+b)))/( np.exp(-1j*a) +  np.exp(-1j*b))

    return val

def test_func2(a,b):

    val = (1 + np.cos((a+b)))/( np.cos(a) +  np.cos(b))

    return val

test_vals = test_func(-2, 5)
test_vals2 = test_func2(-2, 5)

print('test_vals:', test_vals)
print('test_vals2:', test_vals2)

test_vals = test_func(65, 7)
test_vals2 = test_func2(65, 7)
print('test_vals:', test_vals)
print('test_vals2:', test_vals2)

test_vals = test_func(40, 800)
test_vals2 = test_func2(40, 800)
print('test_vals:', test_vals)
print('test_vals2:', test_vals2)

x = y = np.linspace(0, 2*np.pi, 1001)
z = np.array([1/test_func2(i,j) for j in y for i in x])
Z = z.reshape(1001, 1001)

Z[Z<10] = np.nan
#Z[Z<-1] = np.nan

plt.contourf(x, y, Z, 100)

#plt.imshow(Z, interpolation='bilinear')
plt.show()
