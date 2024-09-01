from RIKEN_res_COMSOL_utils import *
import cap_util as cap
import numpy as np
from analytic_coupled_lines_calculator import calc_self_and_coupling_capacitance

## colors
color0= '#ff674b'
color1 = '#c57bff'
color2 = '#ddb2ff'
color3 = '#adb4ff'
color4 = '#d8a4f6'

#initial_lengths = resonator_A.resonator_lengths_from_base_COMSOL_params()

C_val = cap.get_Cc(49e-6)
L_val = cap.get_L(49e-6)
print('C_val:', C_val)


w = 5e-6
s = 7.5e-6
eps_r = 11.7

cs, cm_vals_analytic = calc_self_and_coupling_capacitance([1000e-6], w, s, eps_r)

C_val = cs
C_val = C_val[0]

print('C_val:', C_val)

Z_val = np.sqrt(L_val/C_val)
v_val = 1/(np.sqrt(L_val*C_val))

print('C_val:', C_val)
print('v_val:', v_val)

names = ['A', 'B', 'C', 'D']

### notches V2 pattern

# from RIKEN_V2_double_resonator_pattern_user_params import *

# predicted_notches = []
# predicted_numeric_notches = []

# #names = ['D']

# for name in names:
    
#     readout_params, filter_params = get_res_filter_params(name)
    
#     res_filter_system = RIKEN_coupled_readout_filter_COMSOL(readout_params, filter_params)

#     readout_lengths = res_filter_system.readout.resonator_lengths_from_base_COMSOL_params()
#     #print('readout_lengths:', readout_lengths)
#     total_readout_length = res_filter_system.readout.total_resonator_length_from_base_COMSOL_params()
#     print(name)
#     print('total_readout_length:', total_readout_length)

#     open_lens = res_filter_system.readout.length_open()
#     short_lens = res_filter_system.readout.length_short()
#     coupled_len = res_filter_system.readout.length_coupled()
    
#     print('readout open_lens:', open_lens)
#     print('readout short_lens:', short_lens)
#     print('coupled_len:', coupled_len)

#     open_lens_filter = res_filter_system.filter.length_open()
#     short_lens_filter = res_filter_system.filter.length_short()
#     print('filter open_lens:', open_lens_filter)
#     print('filter short_lens:', short_lens_filter)

#     predicted_notch = res_filter_system.omega_notch(phase_vel = v_val)
#     predicted_numeric_notch = res_filter_system.omega_notch_numeric(Z0 = Z_val, phase_vel = v_val)

#     #print(f'predicted_notch_{name} (GHz):', predicted_notch / (2*np.pi*1e9))
#     #print(f'predicted_numeric_notch_{name} (GHz):', predicted_numeric_notch / (2*np.pi*1e9))

#     J_symb = res_filter_system.J_coupling_symbolic_sol()
#     J_num = res_filter_system.J_coupling_numeric_sol()

#     #print('J_symb:', J_symb / (2*np.pi*1e6))
#     #print('J_num:', J_num / (2*np.pi*1e6))
    
#     omegas = res_filter_system.readout.resonance_omega()

#     #print('res freq (GHz):', omegas/(2*np.pi*1e9))
    
#     predicted_notches.append(predicted_notch)
#     predicted_numeric_notches.append(predicted_numeric_notch)

# sys.exit()

# measured_notches = np.array([8.27, 8.96, 8.69, 8.266]) * 2*np.pi * 1e9 #- 0.1* 2*np.pi * 1e9 ##8.27
# measured_omega_rs = np.array([10264,10690,10479,10050]) * 2*np.pi * 1e6
# measured_omega_ps = np.array([10310,10707,10518,10060]) * 2*np.pi * 1e6

# predicted_numeric_notches = np.array(predicted_numeric_notches)
# predicted_notches = np.array(predicted_notches)

# #plt.plot(predicted_notches/(2*np.pi*1e9), measured_notches/(2*np.pi*1e9), linestyle = 'none', marker = 'o', color = color0, markersize = 8)
# plt.xlim(8, 10)
# plt.ylim(8,10)

# print('notch agreement V2:', (predicted_notches - measured_notches)/measured_notches * 100)

# straight_line = [0,10]
# plt.plot(straight_line, straight_line, linestyle = '--', color = 'k')

### notches V4 pattern

from RIKEN_V4_double_resonator_pattern_user_params import *

predicted_notches = []
predicted_numeric_notches = []

#names = ['D']

for name in names:
    
    readout_params, filter_params = get_res_filter_params(name)
    
    res_filter_system = RIKEN_coupled_readout_filter_COMSOL(readout_params, filter_params)

    readout_lengths = res_filter_system.readout.resonator_lengths_from_base_COMSOL_params()
    #print('readout_lengths:', readout_lengths)
    total_readout_length = res_filter_system.readout.total_resonator_length_from_base_COMSOL_params()
    print(name)
    print('total_readout_length:', total_readout_length)

    open_lens = res_filter_system.readout.length_open()
    short_lens = res_filter_system.readout.length_short()
    coupled_len = res_filter_system.readout.length_coupled()
    
    print('readout open_lens:', open_lens)
    print('readout short_lens:', short_lens)
    print('coupled_len:', coupled_len)

    open_lens_filter = res_filter_system.filter.length_open()
    short_lens_filter = res_filter_system.filter.length_short()
    print('filter open_lens:', open_lens_filter)
    print('filter short_lens:', short_lens_filter)

    predicted_notch = res_filter_system.omega_notch(phase_vel = v_val)
    predicted_numeric_notch = res_filter_system.omega_notch_numeric(Z0 = Z_val, phase_vel = v_val)

    #print(f'predicted_notch_{name} (GHz):', predicted_notch / (2*np.pi*1e9))
    #print(f'predicted_numeric_notch_{name} (GHz):', predicted_numeric_notch / (2*np.pi*1e9))

    # J_symb = res_filter_system.J_coupling_symbolic_sol()
    # J_num = res_filter_system.J_coupling_numeric_sol()

    #print('J_symb:', J_symb / (2*np.pi*1e6))
    #print('J_num:', J_num / (2*np.pi*1e6))
    
    omegas = res_filter_system.readout.resonance_omega()

    #print('res freq (GHz):', omegas/(2*np.pi*1e9))
    
    predicted_notches.append(predicted_notch)
    predicted_numeric_notches.append(predicted_numeric_notch)

#sys.exit()

measured_notches = np.array([8.01, 9.325, 8.775, 8.05]) * 2*np.pi * 1e9 # np.array([8.01, 9.325, 8.775, 8.00]) * 2*np.pi * 1e9 #- 0.1* 2*np.pi * 1e9 ##8.27
measured_omega_rs = np.array([10386,10666,10540,10250]) * 2*np.pi * 1e6
measured_omega_ps = np.array([10407,10710,10566,10232]) * 2*np.pi * 1e6

predicted_numeric_notches = np.array(predicted_numeric_notches)
predicted_notches = np.array(predicted_notches)

plt.plot(predicted_notches/(2*np.pi*1e9), measured_notches/(2*np.pi*1e9), linestyle = 'none', marker = 'o', color = color0, markersize = 8)
plt.xlim(8, 10)
plt.ylim(8,10)

print('notch agreement V4 (to meas):', (predicted_notches - measured_notches)/measured_notches * 100)
print('notch agreement V4 (to sim):', (measured_notches - predicted_notches)/predicted_notches * 100)
avg_notch_agreement_to_sim = np.mean((measured_notches - predicted_notches)/predicted_notches * 100)
print('avg notch agreement (to sim):', avg_notch_agreement_to_sim)

straight_line = [0,10]
plt.plot(straight_line, straight_line, linestyle = '--', color = 'k')

#### plot

ax = plt.gca()
ax.set_aspect('equal', 'box')
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
plt.xlabel(r'Predicted $\omega_n/2\pi$ (GHz)', size = 20)
plt.ylabel(r'Measured $\omega_n/2\pi$ (GHz)', size = 20)
plt.grid(visible=True)
plt.xticks([0,1,2,3,4,5,6,7,8,9,10]) #0,1,2,3,4,5,6,7,8,9,10
plt.yticks([0,1,2,3,4,5,6,7,8,9,10]) # 8,8.5, 9, 9.5, 10
ax.tick_params(axis='both', which='both', width=2.5, labelsize=20)
plt.tight_layout()
plt.savefig('notch_frequeny_agreement.pdf', format = 'pdf')
plt.show()

#########################

# Js
from exact_coupled_transmission_line_eqn_solver import * 

### V2 pattern

# J_preds = []

# measured_notches = np.array([8.27, 8.96, 8.69, 8.266]) * 2*np.pi * 1e9 #- 0.1* 2*np.pi * 1e9 ##8.27
# measured_omega_rs = np.array([10264,10690,10479,10050]) * 2*np.pi * 1e6
# measured_omega_ps = np.array([10310,10707,10518,10060]) * 2*np.pi * 1e6

# for i, name in enumerate(names):
    
#     l_c = 317.5 * 1e-6
#     seps = np.array([5, 5, 5, 5])*1e-6
#     Cm_per_lens = cap.get_Cm(seps)
    
#     J_pred = J_coupling_analytic_by_freqs(measured_omega_rs[i], measured_omega_ps[i], measured_notches[i], l_c, Cm_per_lens[i], phase_vel=v_val, Z0=Z_val, simplified = True)

#     J_preds.append(J_pred)

# J_preds = np.array(J_preds)

# print('J_preds:', J_preds/(2*np.pi*1e6))

# J_measured = np.array([35.1, 30.7, 33.7, 27.9])

# #np.array([2, 2, 3, 1.1])

# #plt.plot(J_preds/(2*np.pi*1e6), J_measured, linestyle = 'none', marker = 'o', color = color0, markersize = 8)
# print('J agreement V2:', (J_preds/(2*np.pi*1e6) - J_measured)/J_measured * 100)

### V4 pattern

J_preds = []
J_preds2 = []

measured_notches = np.array([8.01, 9.325, 8.775, 7.99]) * 2*np.pi * 1e9 #- 0.1* 2*np.pi * 1e9 ##8.27
measured_omega_rs = np.array([10386,10666,10540,10250]) * 2*np.pi * 1e6
measured_omega_ps = np.array([10407,10710,10566,10232]) * 2*np.pi * 1e6

seps = np.array([5.5, 3.8, 4.2, 5.5])*1e-6
Cm_per_lens = cap.get_Cm(seps)
cs, cm_vals_analytic = calc_self_and_coupling_capacitance(seps, w, s, eps_r)
cm_vals_analytic = cm_vals_analytic

c_l = cs + cm_vals_analytic

for i, name in enumerate(names):
    
    l_c = 317.5 * 1e-6
    
    J_pred = J_coupling_analytic_by_freqs(measured_omega_rs[i], measured_omega_ps[i], measured_notches[i], l_c, cm_vals_analytic[i], c_l[i], phase_vel=v_val, simplified = True)

    J_preds.append(J_pred)
    
    avg_omega = (measured_omega_rs[i] + measured_omega_ps[i])/2
    
    k = cm_vals_analytic[i]/c_l[i]

    J2 = 2 * avg_omega/np.pi * k * np.sin(l_c * avg_omega / v_val) * np.cos(np.pi/2 * (avg_omega / measured_notches[i]))
    
    J_preds2.append(J2)

J_preds = np.array(J_preds)
J_preds2 = np.array(J_preds2)

print('J_preds:', J_preds/(2*np.pi*1e6))

J_measured = np.array([39.4, 26.2, 30.9, 36.1])

#np.array([2, 2, 3, 1.1])

plt.plot(J_preds/(2*np.pi*1e6), J_measured, linestyle = 'none', marker = 'o', color = color0, markersize = 8)
print('J agreement V4 (to meas):', (J_preds/(2*np.pi*1e6) - J_measured)/J_measured * 100)
print('J agreement V4 (to sim):', (J_measured - J_preds/(2*np.pi*1e6))/ (J_preds/(2*np.pi*1e6) ) * 100)
avg_J_agreement_to_sim = np.mean((J_measured - J_preds/(2*np.pi*1e6))/ (J_preds/(2*np.pi*1e6) ) * 100)
print('avg J agreement (to sim):', avg_J_agreement_to_sim)
print('J_preds2:', -J_preds2/(2*np.pi*1e6))

#plt.plot(-J_preds2/(2*np.pi*1e6), J_measured, linestyle = 'none', marker = 'o', color = color0, markersize = 8)

### plotting

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
plt.xticks([0,5,10,15,20, 25, 30, 35, 40]) # 0,5,10,15,
plt.yticks([0,5,10,15,20, 25, 30, 35, 40]) #0,5,10,15,
ax.tick_params(axis='both', which='both', width=2.5, labelsize=20)
plt.tight_layout()
plt.savefig('J_coupling_agreement.pdf', format = 'pdf')

plt.show()