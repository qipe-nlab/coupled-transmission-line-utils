from RIKEN_res_COMSOL_utils import *
from RIKEN_V3_double_resonator_pattern_user_params import *

#initial_lengths = resonator_A.resonator_lengths_from_base_COMSOL_params()

names = ['A', 'B', 'C', 'D']

for name in names:
    
    readout_params, filter_params = get_res_filter_params(name)
    
    res_filter_system = RIKEN_coupled_readout_filter_COMSOL(readout_params, filter_params)

    predicted_notch = res_filter_system.omega_notch()
    predicted_numeric_notch = res_filter_system.omega_notch_numeric()

    print(f'predicted_notch_{name} (GHz):', predicted_notch / (2*np.pi*1e9))
    print(f'predicted_numeric_notch_{name} (GHz):', predicted_numeric_notch / (2*np.pi*1e9))

    J_symb = res_filter_system.J_coupling_symbolic_sol()
    #J_num = res_filter_A.J_coupling_numeric_sol()

    print('J_symb:', J_symb / (2*np.pi*1e6))
    #print('J_num:', J_num / (2*np.pi*1e6))

