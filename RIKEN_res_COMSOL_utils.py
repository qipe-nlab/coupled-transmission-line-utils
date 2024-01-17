import numpy as np
import cap_util as cap
from exact_coupled_transmission_line_eqn_solver import * 

class RIKEN_resonator_COMSOL(object):

    def __init__(self, input_dict = None):

        COMSOL_resonator_base_params = \
            {'l_couple_line': 800,
                'l_height': 800,  
                'l_end': 600,
                'l_line': 700, 
                'l_pad_bar': 640,
                'meander_spacing': 100,
                'r': 40,
                'w_line': 20,
                'l_start' : 150,
                'n_meander_curve': 4,
                'start_meander_spacing': 40, ## start meander parameters
                'start_n_meander_curve': 3,
                'start_meander_length': 100,
                'wiggle_offset':20,
                'start_r': 20,
                'first_curve_radius':40,
                'thickness_res_coupler_pad': 37.5

            }

        if input_dict is None:
            input_dict = COMSOL_resonator_base_params

        self.l_couple_line = input_dict['l_couple_line']
        self.l_height = input_dict['l_height']
        self.l_end = input_dict['l_end']
        self.l_line = input_dict['l_line']
        self.l_pad_bar = input_dict['l_pad_bar']
        self.meander_spacing = input_dict['meander_spacing']
        self.r = input_dict['r']
        self.w_line = input_dict['w_line'] # not used
        self.l_start = input_dict['l_start']
        self.n_meander_curve = input_dict['n_meander_curve']
        self.start_meander_spacing = input_dict['start_meander_spacing']
        self.start_n_meander_curve = input_dict['start_n_meander_curve']
        self.start_meander_length = input_dict['start_meander_length']
        self.wiggle_offset = input_dict['wiggle_offset']
        self.start_r = input_dict['start_r']
        self.first_curve_radius = input_dict['first_curve_radius']
        self.thickness_res_coupler_pad = input_dict['thickness_res_coupler_pad']

        self.wiggle_curve_length = self.wiggle_curve_length_func()
        self.start_curve_offset_length = self.start_curve_offset_length_func()
        self.extra_wiggle_length = self.extra_wiggle_len_func()

        return 

    def wiggle_curve_length_func(self):

        wiggle_curve_length_val = np.pi * self.start_meander_spacing/2

        return wiggle_curve_length_val

    def start_curve_offset_length_func(self):
        
        start_curve_offset_length_val = (self.first_curve_radius - self.r) * (2-np.pi/2)
        
        return start_curve_offset_length_val
    
    def extra_wiggle_len_func(self):

        total_wiggle_length = (self.start_meander_length + np.pi * self.start_meander_spacing/2 - self.start_meander_spacing)*self.start_n_meander_curve + np.pi*self.start_r / 2
        wiggle_section_straight_length = self.start_n_meander_curve * self.start_meander_spacing + 2 * self.start_r
        extra_wiggle_len_val = total_wiggle_length - wiggle_section_straight_length

        return extra_wiggle_len_val

    def length_corner_line(self, len_start, len_end, radius, angle = np.pi/2):

        l_start = len_start - radius * np.abs(np.tan(angle/2))

        l_curve = radius * angle 

        l_end = len_end - radius * np.abs(np.tan(angle/2))

        length_corner_line_val = l_start + l_curve + l_end

        return length_corner_line_val

    def length_meander_line(self, spacing, len_each, offset_start, offset_end, n_curves):

        line_length = len_each - spacing
        meander_curve_length = np.pi * spacing/2

        tot_meander_curve_length = meander_curve_length * n_curves
        tot_straight_line_length = line_length * (n_curves-1)

        # print(offset_end)

        start_section_length = 0   #   len_each + offset_start - spacing/2
        end_section_length = line_length/2 + offset_end #   len_each + offset_end - spacing/2

        # print('tot_meander_curve_length:', tot_meander_curve_length)
        # print('tot_straight_line_length:', tot_straight_line_length)
        # print('start_section_length:', start_section_length)
        # print('end_section_length:', end_section_length)
        
        length_meander_line_val = tot_meander_curve_length + tot_straight_line_length + start_section_length + end_section_length

        return length_meander_line_val

    def length_open(self):

        if self.start_n_meander_curve == 0:

            corner_line_1_length = self.length_corner_line(self.l_pad_bar, self.l_start - self.r, self.first_curve_radius, angle = np.pi/2)

            # print('self.l_pad_bar:', self.l_pad_bar)
            # print('corner_line_1_length:', corner_line_1_length)

            start_length = corner_line_1_length - self.thickness_res_coupler_pad
            
            # print('start_length:', start_length)

        else:

            len_start = self.l_pad_bar/2 - self.start_meander_spacing*self.start_n_meander_curve/2 - self.wiggle_offset
            len_end = self.start_r
            radius = self.start_r

            corner_line_1_length = self.length_corner_line(len_start, len_end, radius)

            meander_line_1_length = self.length_meander_line(self.start_meander_spacing, self.start_meander_length, -self.start_r, -self.start_r, self.start_n_meander_curve)

            len_start_2 = self.l_pad_bar/2 - self.start_meander_spacing*self.start_n_meander_curve/2

            corner_line_2_length = self.length_corner_line(len_start_2, len_end, radius)

            len_start_3 = (self.l_pad_bar/2 - self.start_meander_spacing*self.start_n_meander_curve/2)/2 + self.wiggle_offset
            len_end_3 = self.start_r - self.r
            radius_3 = self.first_curve_radius

            corner_line_3_length = self.length_corner_line(len_start_3, len_end_3, radius_3)

            start_length = corner_line_1_length + meander_line_1_length + corner_line_2_length + corner_line_3_length

        mid_length_1 = self.length_corner_line(self.r, self.l_height/2, self.r)
        mid_length_2 = self.length_corner_line(self.l_height/2, self.r, self.r)

        # print('mid_length_1:', mid_length_1)
        # print('mid_length_2:', mid_length_2)

        length_open_val = start_length + mid_length_1 + mid_length_2

        return length_open_val

    def length_short(self):

        start_section_offset_length = self.l_couple_line - self.l_start - self.l_line/2
        # print('start_section_offset_length:', start_section_offset_length)
        
        ## testing: val should be res_1_couple_line_len - res_1_start_len - res_1_line_len/2
        # which is 400 - 190 - 350/2 = 35 # correct
        #start_section_length = 0

        end_section_offset_length = self.l_end - self.l_line/2 - self.extra_wiggle_length +  self.start_curve_offset_length
        
        # print('end_section_offset_length:', end_section_offset_length)
        # print('self.meander_spacing:', self.meander_spacing)

        length_short_val_and_coupled_len = self.length_meander_line(self.meander_spacing, self.l_line, start_section_offset_length, end_section_offset_length, self.n_meander_curve)

        return length_short_val_and_coupled_len

    def length_coupled(self):

        length_coupled_val = self.l_couple_line - self.meander_spacing/2 - self.r

        return length_coupled_val

    def resonator_lengths_from_base_COMSOL_params(self):

        length_open_val = self.length_open()
        length_short_val = self.length_short()
        length_coupled_val = self.length_coupled()

        return length_open_val, length_short_val, length_coupled_val

    def total_resonator_length_from_base_COMSOL_params(self):

        length_open_val, length_short_val, length_coupled_val = self.resonator_lengths_from_base_COMSOL_params()

        total_length = length_open_val + length_short_val + length_coupled_val

        return total_length

    def resonance_omega(self, phase_vel = 1.2e8):

        resonator_length = self.total_resonator_length_from_base_COMSOL_params()*1e-6 # in meters

        lambda_by_4_omega = np.pi*phase_vel/(2*resonator_length)

        return lambda_by_4_omega

    def update_base_params_from_readout_COMSOL_params(self, user_param_dict):

        self.user_param_dict = user_param_dict
        self.resonator_type = 'readout'

        input_dict = self.user_param_dict 

        l_pad_bar_val = input_dict['straight_len'] - input_dict['Q_diam']/2 + input_dict['inner_pad_deformation'] \
        - input_dict['Q_res_coupler_gap'] + input_dict['double_resonator_offset'] - input_dict['double_resonator_spacing'] - 0.95

        l_end_val = input_dict['end_sec_len'] + input_dict['fine_tuning_end'] - input_dict['double_resonator_offset'] 

        self.l_couple_line = input_dict['res_couple_line_len']
        self.l_height = input_dict['height']
        self.l_start = input_dict['start_len']
        self.l_end = l_end_val
        self.l_line = input_dict['l_line_len']
        self.l_pad_bar = l_pad_bar_val 
        self.meander_spacing = input_dict['meander_spacing']
        self.r = input_dict['curvatur_radius']
        self.w_line = input_dict['center_line_width'] # not used
        self.n_meander_curve = input_dict['meander_curves']
        self.start_meander_spacing = input_dict['wiggle_meander_spacing']
        self.start_n_meander_curve = input_dict['wiggle_n']
        self.start_meander_length = input_dict['wiggle_length']
        self.wiggle_offset = input_dict['wiggle_offset']
        self.start_r = input_dict['wiggle_r']
        self.extra_wiggle_length = 0.001122743338823081*1e-6 # hardcoded value in models. Not sure why.
        self.first_curve_radius = input_dict['start_curve_radius']
        self.thickness_res_coupler_pad = input_dict['thickness_res_coupler_pad']
    
        self.wiggle_curve_length = self.wiggle_curve_length_func()
        self.start_curve_offset_length = self.start_curve_offset_length_func()
        self.extra_wiggle_length = self.extra_wiggle_len_func()

        return 
    
    def update_base_params_from_filter_COMSOL_params(self, user_param_dict):
        
        self.user_param_dict = user_param_dict
        self.resonator_type = 'filter'

        input_dict = self.user_param_dict 

        l_pad_bar_val = input_dict['straight_len'] - input_dict['double_resonator_offset'] + input_dict['readout_via_pad_deformation']

        l_end_val = input_dict['end_sec_len'] + input_dict['double_resonator_offset'] -input_dict['readout_via_pad_deformation']

        self.l_couple_line = input_dict['res_couple_line_len']
        self.l_height = input_dict['height']
        self.l_start = input_dict['start_len']
        self.l_end = l_end_val
        self.l_line = input_dict['l_line_len']
        self.l_pad_bar = l_pad_bar_val 
        self.meander_spacing = input_dict['meander_spacing']
        self.r = input_dict['curvatur_radius']
        self.w_line = input_dict['center_line_width'] # not used
        self.n_meander_curve = input_dict['meander_curves']
        self.start_meander_spacing = input_dict['wiggle_meander_spacing']
        self.start_n_meander_curve = input_dict['wiggle_n']
        self.start_meander_length = input_dict['wiggle_length']
        self.wiggle_offset = input_dict['wiggle_offset']
        self.start_r = input_dict['wiggle_r']
        self.extra_wiggle_length = 322.74333882308136e-6 # hardcoded value in models. Not sure why.
        self.first_curve_radius = input_dict['start_curve_radius']
        self.thickness_res_coupler_pad = input_dict['thickness_res_coupler_pad']
    
        self.wiggle_curve_length = self.wiggle_curve_length_func()
        self.start_curve_offset_length = self.start_curve_offset_length_func()
        self.extra_wiggle_length = self.extra_wiggle_len_func()

        return 



## the resonator parameters are now all correct.
    # 

## testing the filter parameters now
    #  173 + 20 + 275
    
resonator_A = RIKEN_resonator_COMSOL()
filter_A = RIKEN_resonator_COMSOL()

#initial_lengths = resonator_A.resonator_lengths_from_base_COMSOL_params()

COMSOL_resonator_user_params_A = \
    {'res_couple_line_len': 400,
        'height': 230,
        'start_len' : 190,
        'end_sec_len': 329,
        'l_line_len': 350,
        'fine_tuning_end': 195, 
        'straight_len': 531.5,
        'Q_diam': 149,
        'inner_pad_deformation': 42.5,
        'Q_res_coupler_gap': 5,
        'double_resonator_offset': 150,
        'double_resonator_spacing': 5.5,
        'meander_spacing': 45,
        'curvatur_radius': 60, # was 40 # maybe 60?
        'center_line_width': 20,
        'meander_curves': 4,
        'wiggle_meander_spacing': 50, ## start meander parameters
        'wiggle_n': 0,
        'wiggle_length': 90,
        'wiggle_offset':60,
        'wiggle_r': 20,
        'start_curve_radius': 110,
        'thickness_res_coupler_pad': 37.5
    }

COMSOL_filter_user_params_A = \
    {'res_couple_line_len': 400,
        'height': 230,
        'start_len' : 190,
        'end_sec_len': 329,
        'l_line_len': 350,
        'fine_tuning_end': 195, 
        'straight_len': 531.5,
        'Q_diam': 149,
        'inner_pad_deformation': 42.5,
        'Q_res_coupler_gap': 5,
        'double_resonator_offset': 150,
        'double_resonator_spacing': 5.5,
        'meander_spacing': 45,
        'curvatur_radius': 60, # was 40 # maybe 60?
        'center_line_width': 20,
        'meander_curves': 4,
        'wiggle_meander_spacing': 35, ## start meander parameters
        'wiggle_n': 0,
        'wiggle_length': 90,
        'wiggle_offset':30,
        'wiggle_r': 20,
        'start_curve_radius': 110,
        'thickness_res_coupler_pad': 32.5,
        'readout_via_pad_deformation': 37.5

    }

resonator_A.update_base_params_from_readout_COMSOL_params(COMSOL_resonator_user_params_A)
filter_A.update_base_params_from_filter_COMSOL_params(COMSOL_filter_user_params_A)

resonator_lengths = resonator_A.resonator_lengths_from_base_COMSOL_params()
filter_lengths = filter_A.resonator_lengths_from_base_COMSOL_params()

print('resonator_lengths:', resonator_lengths)
print('filter_lengths:', filter_lengths)

res_short_len = resonator_lengths[1]
filter_short_len = filter_lengths[1]
res_couple_len = resonator_lengths[2]

filter_frequency = np.pi/2 * 1.2e8/(res_short_len+filter_short_len+res_couple_len) * 1e6
print('filter_frequency (GHz):', filter_frequency/(2*np.pi*1e9))

omega_n = 8.466026082533308*1e9 * 2*np.pi
omega_r = 10.43*1e9 * 2*np.pi
l_c = 317.5*1e-6

predicted_J = J_coupling_analytic_by_freqs(omega_r, omega_r, omega_n, l_c, cap.get_Cm(5.5e-6), phase_vel=3*10**8/2.5, Z0=65)

filter_length = filter_A.total_resonator_length_from_base_COMSOL_params()

C_shunt = 1e-6*filter_length * 1/(1.2e8*65) / 2

print('C_shunt (fF):', C_shunt*1e15)

omega_shift = (100*2*np.pi*1e6/(4*C_shunt*50))**0.5

print('omega_shift (GHZ):', omega_shift/(2*np.pi*1e9))

print('predicted_J (MHz):', predicted_J/(2*np.pi*1e6))

input('.')

test_omega_res_A = resonator_A.resonance_omega()
test_omega_filter_A = filter_A.resonance_omega()

print('test_omega_res_A (GHz):', test_omega_res_A/(2*np.pi*1e9))
print('test_omega_filter_A (GHz):', test_omega_filter_A/(2*np.pi*1e9))