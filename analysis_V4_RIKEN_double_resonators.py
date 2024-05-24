from RIKEN_res_COMSOL_utils import *
from RIKEN_V4_double_resonator_pattern_user_params import *

#initial_lengths = resonator_A.resonator_lengths_from_base_COMSOL_params()

names = ['A', 'B', 'C', 'D']

notches_symbolic = []

for name in names:
    
    readout_params, filter_params = get_res_filter_params(name)
    
    res_filter_system = RIKEN_coupled_readout_filter_COMSOL(readout_params, filter_params)

    readout_lengths = res_filter_system.readout.resonator_lengths_from_base_COMSOL_params()
    #print('readout_lengths:', readout_lengths)
    total_readout_length = res_filter_system.readout.total_resonator_length_from_base_COMSOL_params()
    #print('total_readout_length:', total_readout_length)

    predicted_notch = res_filter_system.omega_notch()
    predicted_numeric_notch = res_filter_system.omega_notch_numeric()

    print(f'predicted_notch_{name} (GHz):', predicted_notch / (2*np.pi*1e9))
    #print(f'predicted_numeric_notch_{name} (GHz):', predicted_numeric_notch / (2*np.pi*1e9))

    J_symb = res_filter_system.J_coupling_symbolic_sol()
    J_num = res_filter_system.J_coupling_numeric_sol()

    #print('J_symb:', J_symb / (2*np.pi*1e6))
    #print('J_num:', J_num / (2*np.pi*1e6))
    notches_symbolic.append(predicted_notch)

##### Predicting from freqs

#sys.exit()

import cap_util as cap

## QA ##

omega_r = 10.385 * 2 * np.pi* 1e9
omega_p = 10.407 * 2 * np.pi* 1e9

omega_n = 8.0 * 2 * np.pi * 1e9
l_c = 317.5 * 1e-6
Cm_per_len = cap.get_Cm(5.5e-6)
print('Cm_per_len:', Cm_per_len/(1e-15))

predicted_J_qub_A = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
print('predicted_J_qub_A (MHz):', predicted_J_qub_A/(2*np.pi*1e6))

## QB ##

omega_r = 10.666 * 2 * np.pi* 1e9
omega_p = 10.712 * 2 * np.pi* 1e9

omega_n = 9.3 * 2 * np.pi * 1e9
l_c = 317.5 * 1e-6
Cm_per_len = cap.get_Cm(3.8e-6)
print('Cm_per_len:', Cm_per_len/(1e-15))

predicted_J_qub_B = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
print('predicted_J_qub_B (MHz):', predicted_J_qub_B/(2*np.pi*1e6))

## QC ##

omega_r = 10.541 * 2 * np.pi* 1e9
omega_p = 10.567 * 2 * np.pi* 1e9

omega_n = 8.8 * 2 * np.pi * 1e9
l_c = 317.5 * 1e-6
Cm_per_len = cap.get_Cm(4.2e-6)
print('Cm_per_len:', Cm_per_len/(1e-15))

predicted_J_qub_C = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
print('predicted_J_qub_C (MHz):', predicted_J_qub_C/(2*np.pi*1e6))

## QD ##

omega_r = 10.247 * 2 * np.pi* 1e9
omega_p = 10.231 * 2 * np.pi* 1e9

omega_n = 8.0 * 2 * np.pi * 1e9
l_c = 317.5 * 1e-6
Cm_per_len = cap.get_Cm(5.5e-6)
print('Cm_per_len:', Cm_per_len/(1e-15))

predicted_J_qub_D = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
print('predicted_J_qub_D (MHz):', predicted_J_qub_D/(2*np.pi*1e6))

measured_notches = np.array([8.0, 9.3, 8.8, 8.0])

notches_symbolic = np.array(notches_symbolic)/(2*np.pi*1e9)

percent_deviation_notches = 100*(notches_symbolic - measured_notches)/measured_notches

print('percent_deviation notches %:', percent_deviation_notches)

line = np.linspace(0, 10, 10) 

plt.plot(line, line, linestyle = '--')
plt.scatter(notches_symbolic, measured_notches)
plt.show()

predicted_Js = np.array([predicted_J_qub_A, predicted_J_qub_B, predicted_J_qub_C, predicted_J_qub_D])/(2*np.pi*1e6)

measured_Js = np.array([39, 27, 31, 36])

percent_deviation_Js = 100*(predicted_Js - measured_Js)/measured_Js

line = np.linspace(0, 40, 10) 

plt.plot(line, line, linestyle = '--')
plt.scatter(predicted_Js, measured_Js)
plt.show()