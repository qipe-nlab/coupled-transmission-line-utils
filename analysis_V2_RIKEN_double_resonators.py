from RIKEN_res_COMSOL_utils import *
from RIKEN_V2_double_resonator_pattern_user_params import *

#initial_lengths = resonator_A.resonator_lengths_from_base_COMSOL_params()

names = ['A', 'B', 'C', 'D']

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

    predicted_notch = res_filter_system.omega_notch()
    predicted_numeric_notch = res_filter_system.omega_notch_numeric()

    print(f'predicted_notch_{name} (GHz):', predicted_notch / (2*np.pi*1e9))
    print(f'predicted_numeric_notch_{name} (GHz):', predicted_numeric_notch / (2*np.pi*1e9))

    J_symb = res_filter_system.J_coupling_symbolic_sol()
    J_num = res_filter_system.J_coupling_numeric_sol()

    print('J_symb:', J_symb / (2*np.pi*1e6))
    print('J_num:', J_num / (2*np.pi*1e6))
    
    omegas = res_filter_system.readout.resonance_omega()

    print('res freq (GHz):', omegas/(2*np.pi*1e9))

measured_notches = np.array([8.27, 8.96, 8.69, 8.27]) * 2*np.pi * 1e9
measured_omega_rs = np.array([10264,10690,10479,10050]) * 2*np.pi * 1e6
measured_omega_ps = np.array([10310,10707,10518,10060]) * 2*np.pi * 1e6

from exact_coupled_transmission_line_eqn_solver import * 

J_preds = []

for i, name in enumerate(names):
    
    l_c = 317.5 * 1e-6
    sep = 5.5*1e-6
    Cm_per_len = cap.get_Cm(sep)
    
    J_pred = J_coupling_analytic_by_freqs(measured_omega_rs[i], measured_omega_ps[i], measured_notches[i], l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)

    J_preds.append(J_pred)

J_preds = np.array(J_preds)

print('J_preds:', J_preds/(2*np.pi*1e6))

##### Predicting from freqs

sys.exit()

import cap_util as cap

## QA ##

omega_r = 10.275 * 2 * np.pi* 1e9
omega_p = 10.275 * 2 * np.pi* 1e9

omega_n = 8.15 * 2 * np.pi * 1e9
l_c = 317.5 * 1e-6
Cm_per_len = cap.get_Cm(5.5e-6)
print('Cm_per_len:', Cm_per_len/(1e-15))

predicted_J_qub_A = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
print('predicted_J_qub_A (MHz):', predicted_J_qub_A/(2*np.pi*1e6))

## QB ##

omega_r = 10.45 * 2 * np.pi* 1e9
omega_p = 10.45 * 2 * np.pi* 1e9

omega_n = 9.1 * 2 * np.pi * 1e9
l_c = 317.5 * 1e-6
Cm_per_len = cap.get_Cm(3.8e-6)
print('Cm_per_len:', Cm_per_len/(1e-15))

predicted_J_qub_B = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
print('predicted_J_qub_B (MHz):', predicted_J_qub_B/(2*np.pi*1e6))

## QC ##

omega_r = 10.39 * 2 * np.pi* 1e9
omega_p = 10.39 * 2 * np.pi* 1e9

omega_n = 8.85 * 2 * np.pi * 1e9
l_c = 317.5 * 1e-6
Cm_per_len = cap.get_Cm(4.2e-6)
print('Cm_per_len:', Cm_per_len/(1e-15))

predicted_J_qub_C = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
print('predicted_J_qub_C (MHz):', predicted_J_qub_C/(2*np.pi*1e6))

## QD ##

omega_r = 10.1 * 2 * np.pi* 1e9
omega_p = 10.1 * 2 * np.pi* 1e9

omega_n = 8.1 * 2 * np.pi * 1e9
l_c = 317.5 * 1e-6
Cm_per_len = cap.get_Cm(5.5e-6)
print('Cm_per_len:', Cm_per_len/(1e-15))

predicted_J_qub_D = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True)
print('predicted_J_qub_D (MHz):', predicted_J_qub_D/(2*np.pi*1e6))