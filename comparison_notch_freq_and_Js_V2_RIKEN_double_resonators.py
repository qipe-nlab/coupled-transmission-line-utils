from RIKEN_res_COMSOL_utils import *
from RIKEN_V2_double_resonator_pattern_user_params import *
import cap_util as cap
import numpy as np

## colors
color0= '#ff674b'
color1 = '#c57bff'
color2 = '#ddb2ff'
color3 = '#adb4ff'
color4 = '#d8a4f6'

#initial_lengths = resonator_A.resonator_lengths_from_base_COMSOL_params()

C_val = cap.get_Cc_plus_Cm(49e-6)
L_val = cap.get_L(49e-6)
Z_val = np.sqrt(L_val/C_val)
v_val = 1/(np.sqrt(L_val*C_val))


names = ['A', 'B', 'C', 'D']
predicted_notches = []
predicted_numeric_notches = []

J_symb_preds = []

names = ['B']

for name in names:
    
    readout_params, filter_params = get_res_filter_params(name)
    
    res_filter_system = RIKEN_coupled_readout_filter_COMSOL(readout_params, filter_params)

    readout_lengths = res_filter_system.readout.resonator_lengths_from_base_COMSOL_params()
    #print('readout_lengths:', readout_lengths)
    total_readout_length = res_filter_system.readout.total_resonator_length_from_base_COMSOL_params()
    print('total_readout_length:', total_readout_length)

    short_lens = res_filter_system.readout.length_short()
    coupled_len = res_filter_system.readout.length_coupled()
    
    print('short_lens:', short_lens)
    print('coupled_len:', coupled_len)

    short_lens_filter = res_filter_system.filter.length_short()
    print('short_lens_filter:', short_lens_filter)

    predicted_notch = res_filter_system.omega_notch(phase_vel = v_val)
    predicted_numeric_notch = res_filter_system.omega_notch_numeric(Z0 = Z_val, phase_vel = v_val)

    print(f'predicted_notch_{name} (GHz):', predicted_notch / (2*np.pi*1e9))
    print(f'predicted_numeric_notch_{name} (GHz):', predicted_numeric_notch / (2*np.pi*1e9))

    J_symb = res_filter_system.J_coupling_symbolic_sol()
    J_num = res_filter_system.J_coupling_numeric_sol()

    print('J_symb:', J_symb / (2*np.pi*1e6))
    print('J_num:', J_num / (2*np.pi*1e6))
    
    omegas = res_filter_system.readout.resonance_omega()

    print('res freq (GHz):', omegas/(2*np.pi*1e9))
    
    predicted_notches.append(predicted_notch)
    predicted_numeric_notches.append(predicted_numeric_notch)

    J_symb_preds.append(J_symb)

measured_notches = np.array([8.27, 8.96, 8.69, 8.266]) * 2*np.pi * 1e9 #- 0.1* 2*np.pi * 1e9 ##8.27
measured_omega_rs = np.array([10264,10690,10479,10050]) * 2*np.pi * 1e6
measured_omega_ps = np.array([10310,10707,10518,10060]) * 2*np.pi * 1e6

J_symb_preds = np.array(J_symb_preds)

###
omega_q = 8.4 * np.pi*1e9

C_q = 70e-15
C_g = 5e-15
C_ext = 19e-15

l_c = coupled_len*1e-6
l_Gf = res_filter_system.readout.length_short()*1e-6
l_Gn = res_filter_system.readout.length_open()*1e-6
l_Rf = res_filter_system.filter.length_short()*1e-6
l_Rn = res_filter_system.filter.length_open()*1e-6

val_g = g_coupling(omega_q, C_q, C_g, l_c, l_Gf, l_Gn, phase_vel=3*10**8/2.5, Z0=65)
val_k = k_filter(C_ext, l_c, l_Rf, l_Rn, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)

print('val_g:', val_g/(2*np.pi*1e6))
print('val_k:', val_k/(2*np.pi*1e6))

omegas = np.linspace(7, 10.5, 100) * 1e9 * 2 * np.pi

d_val = 5e-6
Lm_per_len = cap.get_Lm(d_val) * 1.33
Cm_per_len = cap.get_Cm(d_val) * 1.33


# ## Q3
# omega_r = 10.050*2*np.pi*1e9
# omega_p = 10.060*2*np.pi*1e9
# omega_n = 8.266*2*np.pi*1e9
# J = 30*2*np.pi*1e6

# ## Q2
# omega_r = 10.479*2*np.pi*1e9
# omega_p = 10.501*2*np.pi*1e9
# omega_n = 8.69*2*np.pi*1e9
# J = 33*2*np.pi*1e6

## Q1
omega_r = 10.690*2*np.pi*1e9
omega_p = 10.707*2*np.pi*1e9
omega_n = 8.9*2*np.pi*1e9
J = 38*2*np.pi*1e6

#Q0
# omega_q = 7699 * 2*np.pi * 1e6
# omega_r = 10264 * 2*np.pi * 1e6
# omega_p = 10310 * 2*np.pi * 1e6
# omega_n = 8.3*2*np.pi*1e9
# J = 35 * 2 * np.pi * 1e6

T1_vals = qubit_radiative_decay_equivalent_LE_circuit(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)
T1_vals_no_notch = qubit_radiative_decay_equivalent_LE_circuit_without_notch(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)
T1_from_ham_test = qubit_radiative_decay_from_ham(C_q, C_g, C_ext, omega_r, omega_p, omega_n, J, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)

plt.plot(omegas/(1e9*2*np.pi), T1_vals * 1e3)
plt.plot(omegas/(1e9*2*np.pi), T1_vals_no_notch * 1e3, color ='k')
plt.plot(omegas/(1e9*2*np.pi), T1_from_ham_test * 1e3, color ='k')

np.save('T1_limits_from_circ_Q1', T1_from_ham_test)
np.save('T1_limits_from_circ_Q1_omegas', omegas)

plt.yscale('log')
plt.show()

sys.exit()

## quick_bandwidth_test
alpha = 10

def bandwidth(omega_n, omega_r, omega_p, alpha):
    
    omega_avg = (omega_r + omega_p)/2
    
    val = (1 - (omega_n/omega_avg)**2)/alpha**0.5 * omega_n
    
    return val

B_val = bandwidth(measured_notches[0], measured_omega_rs[0], measured_omega_ps[0], alpha)

print('B_val (MHz):', B_val/(2*np.pi*1e6))

#sys.exit()

from exact_coupled_transmission_line_eqn_solver import * 

J_preds = []

for i, name in enumerate(names):
    
    l_c = 337.5 * 1e-6
    sep = 5*1e-6
    Cm_per_len = cap.get_Cm(sep)
    
    J_pred = J_coupling_analytic_by_freqs(measured_omega_rs[i], measured_omega_ps[i], measured_notches[i], l_c, Cm_per_len, phase_vel=v_val, Z0=Z_val, simplified = True)

    J_preds.append(J_pred)

J_preds = np.array(J_preds)

print('J_preds:', J_preds/(2*np.pi*1e6))

#### testing error on Js

# ##### Predicting from freqs

# #sys.exit()

# import cap_util as cap

# ## QA ##

# omega_r = 10.275 * 2 * np.pi* 1e9
# omega_p = 10.275 * 2 * np.pi* 1e9

# omega_n = 8.15 * 2 * np.pi * 1e9
# l_c = 317.5 * 1e-6
# Cm_per_len = cap.get_Cm(5.5e-6)
# print('Cm_per_len:', Cm_per_len/(1e-15))

# predicted_J_qub_A = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
# print('predicted_J_qub_A (MHz):', predicted_J_qub_A/(2*np.pi*1e6))

# ## QB ##

# omega_r = 10.45 * 2 * np.pi* 1e9
# omega_p = 10.45 * 2 * np.pi* 1e9

# omega_n = 9.1 * 2 * np.pi * 1e9
# l_c = 317.5 * 1e-6
# Cm_per_len = cap.get_Cm(3.8e-6)
# print('Cm_per_len:', Cm_per_len/(1e-15))

# predicted_J_qub_B = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
# print('predicted_J_qub_B (MHz):', predicted_J_qub_B/(2*np.pi*1e6))

# ## QC ##

# omega_r = 10.39 * 2 * np.pi* 1e9
# omega_p = 10.39 * 2 * np.pi* 1e9

# omega_n = 8.85 * 2 * np.pi * 1e9
# l_c = 317.5 * 1e-6
# Cm_per_len = cap.get_Cm(4.2e-6)
# print('Cm_per_len:', Cm_per_len/(1e-15))

# predicted_J_qub_C = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
# print('predicted_J_qub_C (MHz):', predicted_J_qub_C/(2*np.pi*1e6))

# ## QD ##

# omega_r = 10.1 * 2 * np.pi* 1e9
# omega_p = 10.1 * 2 * np.pi* 1e9

# omega_n = 8.1 * 2 * np.pi * 1e9
# l_c = 317.5 * 1e-6
# Cm_per_len = cap.get_Cm(5.5e-6)
# print('Cm_per_len:', Cm_per_len/(1e-15))

# predicted_J_qub_D = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
# print('predicted_J_qub_D (MHz):', predicted_J_qub_D/(2*np.pi*1e6))

###### 

J_measured = np.array([35.1, 30.7, 33.7, 27.9])

test_Js = np.array([[38.6, 32., 35.7, 28.1],[35.7,34.6,37.1,30.5], [34.4,36.0,32.7,27.9],[42.7,29.5,33.9,30.4], [37.1,33.8,33,26.7], [40.6,30.8,35.0,29.3]])

#test_Js_measured = np.mean(test_Js, axis = 0)
#test_Js_errors  = np.std(test_Js, axis = 0)

#print('test_Js_measured:', test_Js_measured)
#print('test_Js_errors:', test_Js_errors)

J_errors = np.array([0.1, 0.1, 0.1, 0.1])

#np.array([2, 2, 3, 1.1])

#plt.figure(figsize=(6, 6))
plt.errorbar(J_preds/(2*np.pi*1e6), J_measured, yerr = J_errors, linestyle = 'none', marker = 'o', color = color0, markersize = 8)
#plt.errorbar(J_preds/(2*np.pi*1e6), test_Js_measured, yerr = test_Js_errors, linestyle = 'none', marker = 'o', color = 'k', markersize = 8)
#plt.errorbar(J_symb_preds/(2*np.pi*1e6), J_measured, yerr = 3, xerr = 1, linestyle = 'none', marker = 'o', color = 'purple', markersize = 8)

print('J agreement:', (J_preds/(2*np.pi*1e6) - J_measured)/J_measured * 100)

plt.xlim(20, 40)
plt.ylim(20,40)

straight_line = [0,50]

plt.plot(straight_line, straight_line, linestyle = '--', color = 'k')

ax = plt.gca()
ax.set_aspect('equal', 'box')
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
plt.xlabel(r'Predicted $J/2\pi$ (MHz)', size = 20)
plt.ylabel(r'Measured $J/2\pi$ (MHz)', size = 20)
plt.grid(visible=True)
plt.xticks([20, 25, 30, 35, 40]) # 0,5,10,15,
plt.yticks([20, 25, 30, 35, 40]) #0,5,10,15,
ax.tick_params(axis='both', which='both', width=2.5, labelsize=20)

plt.tight_layout()
plt.show()

#### notches

predicted_numeric_notches = np.array(predicted_numeric_notches)
predicted_notches = np.array(predicted_notches)
predicted_notches_x_err = np.array([0,0,0,0])
predicted_notches_y_err = np.array([100, 20, 50, 20])*1e-3

#plt.errorbar(predicted_numeric_notches/(2*np.pi*1e9), measured_notches/(2*np.pi*1e9),  linestyle = 'none', marker = 'o', color = 'blue')
plt.errorbar(predicted_notches/(2*np.pi*1e9), measured_notches/(2*np.pi*1e9), yerr =predicted_notches_y_err , xerr = predicted_notches_x_err/(2*np.pi*1e9)*0.02, linestyle = 'none', marker = 'o', color = color0, markersize = 8)
plt.xlim(8, 10)
plt.ylim(8,10)

print('notch agreement:', (predicted_notches - measured_notches)/measured_notches * 100)

straight_line = [0,10]
plt.plot(straight_line, straight_line, linestyle = '--', color = 'k')

ax = plt.gca()
ax.set_aspect('equal', 'box')
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
plt.xlabel(r'Predicted $\omega_n/2\pi$ (GHz)', size = 20)
plt.ylabel(r'Measured $\omega_n/2\pi$ (GHz)', size = 20)
plt.grid(visible=True)
plt.xticks([8,8.5, 9, 9.5, 10])
plt.yticks([8,8.5, 9, 9.5, 10])
ax.tick_params(axis='both', which='both', width=2.5, labelsize=20)

plt.tight_layout()
plt.show()