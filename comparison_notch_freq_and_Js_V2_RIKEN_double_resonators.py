from RIKEN_res_COMSOL_utils import *
from RIKEN_V2_double_resonator_pattern_user_params import *
import cap_util as cap

#initial_lengths = resonator_A.resonator_lengths_from_base_COMSOL_params()

C_val = cap.get_Cc_plus_Cm(49e-6)
L_val = cap.get_L(49e-6)
Z_val = np.sqrt(L_val/C_val)
v_val = 1/(np.sqrt(L_val*C_val))


names = ['A', 'B', 'C', 'D']
predicted_notches = []
predicted_numeric_notches = []

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

measured_notches = np.array([8.27, 8.96, 8.69, 8.266]) * 2*np.pi * 1e9 ##8.27
measured_omega_rs = np.array([10264,10690,10479,10050]) * 2*np.pi * 1e6
measured_omega_ps = np.array([10310,10707,10518,10060]) * 2*np.pi * 1e6

## quick_bandwidth_test
alpha = 10

def bandwidth(omega_n, omega_r, omega_p, alpha):
    
    omega_avg = (omega_r + omega_p)/2
    
    val = (1 - (omega_n/omega_avg)**2)/alpha**0.5 * omega_n
    
    return val

B_val = bandwidth(measured_notches[0], measured_omega_rs[0], measured_omega_ps[0], alpha)

print('B_val (MHz):', B_val/(2*np.pi*1e6))

sys.exit()

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

J_measured = np.array([33.1,26.2,28.3,21.4])
#J_errors = 

#plt.figure(figsize=(6, 6))

plt.errorbar(J_preds/(2*np.pi*1e6), J_measured, yerr = 3, xerr = 1, linestyle = 'none', marker = 'o', color = 'purple', markersize = 8)
plt.xlim(20, 35)
plt.ylim(20,35)

straight_line = [0,35]
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
plt.xticks([20,25,30,35]) # 0,5,10,15,
plt.yticks([20,25,30,35]) #0,5,10,15,
ax.tick_params(axis='both', which='both', width=2.5, labelsize=20)

plt.tight_layout()
plt.show()

#### notches

predicted_numeric_notches = np.array(predicted_numeric_notches)
predicted_notches = np.array(predicted_notches)

#plt.errorbar(predicted_numeric_notches/(2*np.pi*1e9), measured_notches/(2*np.pi*1e9),  linestyle = 'none', marker = 'o', color = 'blue')
plt.errorbar(predicted_notches/(2*np.pi*1e9), measured_notches/(2*np.pi*1e9), yerr = 0.05, xerr = 0.1, linestyle = 'none', marker = 'o', color = 'purple', markersize = 8)
plt.xlim(8, 10)
plt.ylim(8,10)

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