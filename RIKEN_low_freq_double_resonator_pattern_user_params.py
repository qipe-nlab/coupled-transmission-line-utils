
COMSOL_readout_filter_shared_params = \
    {   'res_couple_line_len': 750, # 400
        'height': 230,
        'start_len' : 375,    
        'straight_len': 531.5,
        'l_line_len': 700,
        'Q_res_coupler_gap': 5,
        'double_resonator_offset': 125,
        'meander_spacing': 45,
        'curvatur_radius': 60,
        'center_line_width': 20,
        'readout_via_pad_deformation': 37.5 
    }

COMSOL_readout_shared_params = \
    {   'wiggle_offset':60,
        'wiggle_length': 90,
        'wiggle_n':0,
        'wiggle_r': 20,
        'wiggle_meander_spacing': 50,
     }

COMSOL_filter_shared_params = \
    {   'wiggle_offset':30,
        'wiggle_r': 20,
        'wiggle_meander_spacing': 35,
        'thickness_res_coupler_pad': 32.5,
        'readout_via_pad_deformation': 37.5
     }

### unique parameters begin here

COMSOL_readout_user_params_A = \
    {   'end_sec_len': 100,
        'fine_tuning_end': 178, 
        'Q_diam': 270,
        'inner_pad_deformation': 35,
        'double_resonator_spacing': 8,
        'start_curve_radius': 110,
        'meander_curves': 5,
    }

COMSOL_filter_user_params_A = \
    {   'wiggle_n': 0,
        'wiggle_length': 90,
        'meander_curves': 5,
    }

COMSOL_readout_user_params_B = \
    {   'end_sec_len': 400,
        'fine_tuning_end': 150, 
        'Q_diam': 240,
        'inner_pad_deformation': 17.5,
        'double_resonator_spacing': 6, #3.8
        'start_curve_radius': 60,
        'meander_curves': 4,
    }

COMSOL_filter_user_params_B = \
    {   'wiggle_n': 3,
        'wiggle_length': 82.5,
        'meander_curves': 4,
    }

COMSOL_readout_user_params_C = \
    {   'end_sec_len': 450,
        'fine_tuning_end': 150, 
        'Q_diam': 250,
        'inner_pad_deformation': 17.5,
        'double_resonator_spacing': 6.5, # 4.2
        'start_curve_radius': 60,
        'meander_curves': 4,
    }

COMSOL_filter_user_params_C = \
    {   'wiggle_n': 3,
        'wiggle_length': 85,
        'meander_curves': 4,
    }

COMSOL_readout_user_params_D = \
    {   'end_sec_len': 130,
        'fine_tuning_end': 220, 
        'Q_diam': 280,
        'inner_pad_deformation': 35,
        'double_resonator_spacing': 7.5,
        'start_curve_radius': 130,
        'meander_curves': 5,
    }

COMSOL_filter_user_params_D = \
    {   'wiggle_n': 0,
        'wiggle_length': 90,
        'meander_curves': 5,
    }

def generate_readout_param_dict(unit):

    shared_res_dict = COMSOL_readout_filter_shared_params | COMSOL_readout_shared_params

    if unit == 'A':
        param_dict = COMSOL_readout_user_params_A
    if unit == 'B':
        param_dict = COMSOL_readout_user_params_B
    if unit == 'C':
        param_dict = COMSOL_readout_user_params_C
    if unit == 'D':
        param_dict = COMSOL_readout_user_params_D

    total_dict = shared_res_dict | param_dict

    return total_dict

def generate_filter_param_dict(unit):

    shared_res_dict = COMSOL_readout_filter_shared_params | COMSOL_filter_shared_params

    if unit == 'A':
        param_dict_readout = COMSOL_readout_user_params_A
        param_dict_filter = COMSOL_filter_user_params_A
    if unit == 'B':
        param_dict_readout = COMSOL_readout_user_params_B
        param_dict_filter = COMSOL_filter_user_params_B

    if unit == 'C':
        param_dict_readout = COMSOL_readout_user_params_C
        param_dict_filter = COMSOL_filter_user_params_C
                             
    if unit == 'D':
        param_dict_readout = COMSOL_readout_user_params_D
        param_dict_filter = COMSOL_filter_user_params_D

    param_dict_filter = param_dict_readout| param_dict_filter 

    total_dict = shared_res_dict | param_dict_filter

    return total_dict

readout_A_params = generate_readout_param_dict('A')
readout_B_params = generate_readout_param_dict('B')
readout_C_params = generate_readout_param_dict('C')
readout_D_params = generate_readout_param_dict('D')

filter_A_params = generate_filter_param_dict('A')
filter_B_params = generate_filter_param_dict('B')
filter_C_params = generate_filter_param_dict('C')
filter_D_params = generate_filter_param_dict('D')

def get_res_filter_params(unit):

    if unit == 'A':
        readout_params = readout_A_params
        filter_params = filter_A_params
    if unit == 'B':
        readout_params = readout_B_params
        filter_params = filter_B_params
    if unit == 'C':
        readout_params = readout_C_params
        filter_params = filter_C_params
    if unit == 'D':
        readout_params = readout_D_params
        filter_params = filter_D_params

    return readout_params, filter_params