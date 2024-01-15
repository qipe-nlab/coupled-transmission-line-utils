from cmath import phase
import re
import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings

#######################
### Basic functions ###
#######################

def Zcap(Cval, omega):

    if Cval == 0:
        val = 1j*1e30
    else:
        val = -1j/(Cval*omega)

    return val

def Zind(Lval, omega):
    return 1j*(Lval*omega)

def Zopen(Z0, phi):
    return Z0/(1j*np.tan(phi))

def Zshort(Z0, phi):
    return 1j*Z0*np.tan(phi)

def Zinput(Z0, Zload, phi):
    return Z0*(Zload+1j*Z0*np.tan(phi))/(Z0+1j*Zload*np.tan(phi))

def Zpara(Z1, Z2):
    return 1/(1/Z1+1/Z2)

def Zres(Cr, Lr, omega):
    return Zpara(Zcap(Cr, omega),Zind(Lr, omega))

def Zchar(Cr, Lr):
    return (Lr/Cr)**0.5

def omega_r(Cr, Lr):
    return 1/(Lr*Cr)**0.5

########################################
### Distributed Element Calculations ###
###     General Numerical Solver     ###
########################################

# General exact solver for n coupled lines
# Only assumption is that lines are lossless
# See: Analysis of Multiconductor Transmission Lines, 2E, Chapter 7

# defining matricies

def Z_mat(L_matrix, omega):

    mat = 1j*omega*L_matrix

    return mat

def Y_mat(C_matrix, omega):

    mat = 1j*omega*C_matrix

    return mat 

def C_mutual_to_Maxwell(mutual_C_matrix):

    diags = np.sum(mutual_C_matrix, axis = 1)

    mat = - mutual_C_matrix.copy()

    np.fill_diagonal(mat, diags)

    return mat 

def U_mat_and_eig(C_matrix):

    # for diagonalizing C

    eig, mat = np.linalg.eigh(C_matrix)

    return mat, eig

def S_mat_and_eig(U_mat, U_eig, L_mat):

    # for diagonalizing L

    sqrt_U_eig = np.emath.sqrt(U_eig)

    sqrt_U_eig_mat = np.diag(sqrt_U_eig)

    M_mat = sqrt_U_eig_mat @ (U_mat.T) @ L_mat @ U_mat @ sqrt_U_eig_mat

    eig, mat = np.linalg.eigh(M_mat)

    return mat, eig

def T_mat(U_mat, S_mat, U_eig):

    # for diagonalizing CL matrix

    sqrt_U_eig = np.emath.sqrt(U_eig)

    mat = U_mat@np.diag(sqrt_U_eig)@S_mat

    return mat

def T_v_mat(U_mat, S_mat, U_eig):

    # for diagonalizing CL matrix

    sqrt_U_eig = np.emath.sqrt(U_eig)

    mat = U_mat@np.diag(1/sqrt_U_eig)@S_mat

    return mat

def T_inv_mat(U_mat, S_mat, U_eig):

    mat = (T_v_mat(U_mat, S_mat, U_eig)).T

    return mat 

def propagation_eig(S_eig, omega):

    # diagonalized YZ matrix elements

    # print('omega:', omega)
    # print('S_eig:', S_eig)

    eig = 1j*omega*np.emath.sqrt(S_eig)

    #  print('eig:', eig)

    return eig

def Z_char_mat(U_mat, S_mat, U_eig, S_eig):

    sqrt_U_eig = np.emath.sqrt(U_eig)

    sqrt_S_eig = np.emath.sqrt(S_eig)

    mat = U_mat@np.diag(1/sqrt_U_eig)@S_mat@np.diag(sqrt_S_eig)@(S_mat.T)@np.diag(1/sqrt_U_eig)@(U_mat.T)

    return mat

def prop_exp_mats(prop_eigs, pos):

    forward_prop_mat = np.diag(np.exp(-prop_eigs*pos))

    backward_prop_mat = np.diag(np.exp(prop_eigs*pos))

    return forward_prop_mat, backward_prop_mat

def Z_near_bound_mat(Z_near_list):

    # boundary conditions on near side

    mat = np.diag(Z_near_list)

    return mat

def Z_far_bound_mat(Z_far_list):

    # boundary conditions on far side

    mat = np.array(Z_far_list)

    return mat

def defining_mat(Zchar, Z_nb, Z_fb, T, prop_eigs, len):

    # defining matrix for determining current solutions. 
    # Eq. 7.90 of above reference.
    
    # def: defining_mat@(uc_I_vec) = V_vec

    # nb: near bound
    # fb: far bound

    # print('Zchar, Z_nb:', Zchar, Z_nb)
    # print('Z_fb:', Z_fb)
    # print('T:', T)
    # print('prop_eigs:', prop_eigs)
    # print('len:', len)

    fp_exp_end, bp_exp_end = prop_exp_mats(prop_eigs, len)

    mat11 = (Zchar + Z_nb)@T
    mat12 = (Zchar - Z_nb)@T
    mat21 = (Zchar - Z_fb)@T@fp_exp_end
    mat22 = (Zchar + Z_fb)@T@bp_exp_end

    mat = np.block([[mat11, mat12],[mat21, mat22]])

    return mat

def uc_I_sols(defining_matrix, Vs_nb, Vs_fb):

    # solve for the diagonalized directional currents along the lines given the source voltages
    # achieved by inverting the defining matrix (eq. 7.90)

    # print('defining_matrix:', defining_matrix)
    # print('type(defining_matrix):', type(defining_matrix))

    defining_mat_inv = np.linalg.inv(defining_matrix)

    Vs = np.concatenate((Vs_nb, Vs_fb))

    Is = defining_mat_inv@Vs

    Is_split = np.array_split(Is, 2)

    forward_uc_I_sols = Is_split[0]
    backward_uc_I_sols = Is_split[1]

    return forward_uc_I_sols, backward_uc_I_sols 

def V_I_sols(Z_char_mat_sol, T_mat_sol, forward_uc_I_sols, backward_uc_I_sols, prop_eigs, pos):

    # solve for the voltage and current solutions along the lines in the original (undiagonalized) frame
    # achieved by subbing in the directional current solutions into eqs. 7.86
    # uc : uncoupled

    forward_prop_mat, backward_prop_mat = prop_exp_mats(prop_eigs, pos)

    sum_current_vec = forward_prop_mat@forward_uc_I_sols + backward_prop_mat@backward_uc_I_sols
    
    diff_current_vec = forward_prop_mat@forward_uc_I_sols - backward_prop_mat@backward_uc_I_sols

    V_sols = Z_char_mat_sol@T_mat_sol@sum_current_vec

    I_sols = T_mat_sol@diff_current_vec

    return V_sols, I_sols

# defining matricies for special case of symmetric 2 line problem

def T_mat_sym_three():

    # for diagonalizing CL matrix

    mat = 1/(2)**0.5 * np.array([[1, 1],[1, -1]])

    return mat

def Z_char_mat_sym_three(Z_mat, Y_mat):

    Z_char_plus = np.emath.sqrt((Z_mat[0,0] + Z_mat[0,1])/(Y_mat[0,0] + Y_mat[0,1]))
    Z_char_minus =  np.emath.sqrt((Z_mat[0,0] - Z_mat[0,1])/(Y_mat[0,0] - Y_mat[0,1]))

    mat = 0.5 * np.array([[Z_char_plus + Z_char_minus, Z_char_plus - Z_char_minus],[Z_char_plus - Z_char_minus, Z_char_plus + Z_char_minus]])

    return mat

def propagation_eig_sym_three(Z_mat, Y_mat):

    # diagonalized YZ matrix elements

    eig1 = np.emath.sqrt((Z_mat[0,0] + Z_mat[0,1])*(Y_mat[0,0] + Y_mat[0,1]))
    
    eig2 = np.emath.sqrt((Z_mat[0,0] - Z_mat[0,1])*(Y_mat[0,0] - Y_mat[0,1]))

    eig = np.array([eig1, eig2])

    return eig

# basic util functions

def Z_short_tl(Z0, phase_vel, L, omega):

    # impedance of a shorted TL section 

    tau = L/phase_vel

    val = 1j*Z0*np.tan(tau * omega)

    return val 

def Z_open_tl(Z0, phase_vel, L, omega):

    # impedance of an open TL section 

    tau = L/phase_vel

    val = -1j*Z0/np.tan(tau * omega)

    return val

def Z_input_tl(Z0, Zload, phase_vel, L, omega ):

    phi = L * omega/phase_vel

    return Z0*(Zload+1j*Z0*np.tan(phi))/(Z0+1j*Zload*np.tan(phi))

def voltage_at_source_location_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, Z0, phase_vel, omega, solved_Z_input_coupled_system = None):

    # returns the exact ratio of the input current to the generator line to the voltage at the start of the coupled line section.

    # if a solved value for the input impedance of the coupled system is provided, then use it directly, otherwise solve it here first

    if solved_Z_input_coupled_system is None:
        Zinput_total = Z_input_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega, phase_vel=phase_vel, Z0=Z0, Zline = 50, C_env = 0)

    else:
        Z_input_coupled_system = solved_Z_input_coupled_system

        Zinput_total = Z_input_tl(Z0, Z_input_coupled_system, phase_vel, l_Gn, omega)

    I_in = 1

    v_in = Zinput_total * I_in

    tau_n = l_Gn / phase_vel

    v_out, I_out = transmission_line_voltage_current_out(Z0, tau_n*omega, v_in, I_in)

    val = v_out

    return val

# Solvers - for general case and for special case of 2 symmetric coupled lines

def V_I_solutions_general(C_mat, L_mat, Z_nb, Z_fb, Vs_nb, Vs_fb, len, omega):

    U_mat_sol, U_eig_sol = U_mat_and_eig(C_mat)
    
    S_mat_sol, S_eig_sol = S_mat_and_eig(U_mat_sol, U_eig_sol, L_mat)

    T_mat_sol = T_mat(U_mat_sol, S_mat_sol, U_eig_sol)

    propagation_eig_sol = propagation_eig(S_eig_sol, omega)

    Z_char_mat_sol = Z_char_mat(U_mat_sol, S_mat_sol, U_eig_sol, S_eig_sol) # L_mat

    defining_mat_sol = defining_mat(Z_char_mat_sol, Z_nb, Z_fb, T_mat_sol, propagation_eig_sol, len)
    
    #print('a')

    forward_uc_I_sols, backward_uc_I_sols  = uc_I_sols(defining_mat_sol, Vs_nb, Vs_fb)

    V_sols_n, I_sols_n = V_I_sols(Z_char_mat_sol, T_mat_sol, forward_uc_I_sols, backward_uc_I_sols, propagation_eig_sol, 0)

    V_sols_f, I_sols_f = V_I_sols(Z_char_mat_sol, T_mat_sol, forward_uc_I_sols, backward_uc_I_sols, propagation_eig_sol, len)

    return V_sols_n, V_sols_f, I_sols_n, I_sols_f

def V_I_solutions_sym_three(C_mat, L_mat, Z_nb, Z_fb, Vs_nb, Vs_fb, len, omega):

    Z_matrix_sol = Z_mat(L_mat, omega)
    Y_matrix_sol = Y_mat(C_mat, omega)

    T_mat_sol = T_mat_sym_three()

    propagation_eig_sol = propagation_eig_sym_three(Z_matrix_sol, Y_matrix_sol)

    Z_char_mat_sol = Z_char_mat_sym_three(Z_matrix_sol, Y_matrix_sol)

    defining_mat_sol = defining_mat(Z_char_mat_sol, Z_nb, Z_fb, T_mat_sol, propagation_eig_sol, len)

    #print('b')

    forward_uc_I_sols, backward_uc_I_sols  = uc_I_sols(defining_mat_sol, Vs_nb, Vs_fb)

    V_sols_n, I_sols_n = V_I_sols(Z_char_mat_sol, T_mat_sol, forward_uc_I_sols, backward_uc_I_sols, propagation_eig_sol, 0)

    V_sols_f, I_sols_f = V_I_sols(Z_char_mat_sol, T_mat_sol, forward_uc_I_sols, backward_uc_I_sols, propagation_eig_sol, len)

    return V_sols_n, V_sols_f, I_sols_n, I_sols_f

########################################
### Distributed Element Calculations ###
### 2 coupled transmission lines, WC ###
########################################

### For 2 coupled transmission lines under the assumption of weak coupling
### See: Solution of the Transmission-Line Equations Under the Weak-Coupling Assumption - C.R. Paul

def transmission_line_voltage_current_out(Z0, electric_length, voltage_in, current_in):

    voltage_out = voltage_in*np.cos(electric_length)-1j*current_in*Z0*np.sin(electric_length)

    current_out = -1j*voltage_in/Z0*np.sin(electric_length)+current_in*np.cos(electric_length)

    return voltage_out, current_out

def A_val(k_plus, k_minus, tau, gamma_Gf, omega):

    return k_plus/(2*tau) * (1-np.exp(-1j*omega*2*tau)) - gamma_Gf * k_minus * np.exp(-1j*omega*2*tau) * 1j*omega

def B_val(k_plus, k_minus, tau, gamma_Gf, omega):

    return gamma_Gf * k_plus/(2*tau) * (1-np.exp(-1j*omega*2*tau)) * np.exp(-1j*omega*tau) - k_minus * 1j*omega * np.exp(-1j*omega*tau)

def reflection_ceofficient(ZL, Z0):

    return (ZL - Z0)/(ZL + Z0)

def F_val(gamma1, gamma2, tau, omega):

    return 1/(1-gamma1*gamma2*np.exp(-1j*omega*2*tau))

def voltage_transmission_coupled_lines_debug(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega):

    ## induced voltage at near side of reciever line given a 1 volt excitation at the generator near side.

    Lm = Lm_per_len * l_c
    Cm = Cm_per_len * l_c

    #### testing using exact expressions:
    C_tl, L_tl = transmission_line_C_and_L(phase_vel, Z0)
    Z0_c = Zchar(C_tl + Cm_per_len, L_tl)
    phase_vel_c = omega_r(C_tl + Cm_per_len, L_tl)

    K_plus = Lm/Z0_c + Cm*Z0_c

    K_minus = Lm/Z0_c - Cm*Z0_c

    tau_c = l_c / phase_vel_c

    Z_Rn = Z_open_tl(Z0, phase_vel, l_Rn, omega)
    Z_Rf = Z_short_tl(Z0, phase_vel, l_Rf, omega)
    Z_Gf = Z_short_tl(Z0, phase_vel, l_Gf, omega)

    l_G = l_c +l_Gf + l_Gn

    tau_G = l_G / phase_vel
    tau_Gn = l_Gn / phase_vel

    ## test changes
    Z_Gn = 0 # 0 # 1j * Z0 * (np.tan(tau_G*omega)*np.cos(tau_Gn*omega) - np.sin(tau_Gn*omega))  / (np.tan(tau_G*omega)*np.sin(tau_Gn*omega) + np.cos(tau_Gn*omega)) 

    gamma_Rn = reflection_ceofficient(Z_Rn, Z0)

    gamma_Rf = reflection_ceofficient(Z_Rf, Z0)

    gamma_Gn = reflection_ceofficient(Z_Gn, Z0)

    gamma_Gf = reflection_ceofficient(Z_Gf, Z0)

    F_G = F_val(gamma_Gn, gamma_Gf, tau_c, omega)

    F_R = F_val(gamma_Rn, gamma_Rf, tau_c, omega)

    A_value = A_val(K_plus, K_minus, tau_c, gamma_Gf, omega)

    B_value = B_val(K_plus, K_minus, tau_c, gamma_Gf, omega)

    #val = (Z_Rn/(Z_Rn + Z0_c)) * (A_value + B_value * gamma_Rf * np.exp(-1j*omega*tau_c)) * F_G * F_R * (Z0_c/(Z0_c + Z_Gn))
    
    ## just return the part that goes to 0 at the notch
    val = (A_value + B_value * gamma_Rf * np.exp(1j*omega*tau_c))

    #val = (A_value + B_value * gamma_Rf * np.exp(1j*omega*tau_c)) *  F_G * F_R

    return val

def voltage_transmission_coupled_lines(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega):

    ## induced voltage at near side of reciever line given a 1 volt excitation at the generator near side.

    Lm = Lm_per_len * l_c
    Cm = Cm_per_len * l_c

    #### testing using exact expressions:
    C_tl, L_tl = transmission_line_C_and_L(phase_vel, Z0)
    Z0_c = Zchar(C_tl + Cm_per_len, L_tl)
    phase_vel_c = omega_r(C_tl + Cm_per_len, L_tl)

    K_plus = Lm/Z0_c + Cm*Z0_c

    K_minus = Lm/Z0_c - Cm*Z0_c

    tau_c = l_c / phase_vel_c

    Z_Rn = Z_open_tl(Z0, phase_vel, l_Rn, omega)
    Z_Rf = Z_short_tl(Z0, phase_vel, l_Rf, omega)
    Z_Gf = Z_short_tl(Z0, phase_vel, l_Gf, omega)

    l_G = l_c +l_Gf + l_Gn

    tau_G = l_G / phase_vel
    tau_Gn = l_Gn / phase_vel

    # test changes
    Z_Gn = 0 #0 # 1j * Z0 * (np.tan(tau_G*omega)*np.cos(tau_Gn*omega) - np.sin(tau_Gn*omega))  / (np.tan(tau_G*omega)*np.sin(tau_Gn*omega) + np.cos(tau_Gn*omega)) 

    gamma_Rn = reflection_ceofficient(Z_Rn, Z0)

    gamma_Rf = reflection_ceofficient(Z_Rf, Z0)

    gamma_Gn = reflection_ceofficient(Z_Gn, Z0)

    gamma_Gf = reflection_ceofficient(Z_Gf, Z0)

    F_G = F_val(gamma_Gn, gamma_Gf, tau_c, omega)

    F_R = F_val(gamma_Rn, gamma_Rf, tau_c, omega)

    A_value = A_val(K_plus, K_minus, tau_c, gamma_Gf, omega)

    B_value = B_val(K_plus, K_minus, tau_c, gamma_Gf, omega)

    val = (Z_Rn/(Z_Rn + Z0_c)) * (A_value + B_value * gamma_Rf * np.exp(-1j*omega*tau_c)) * F_G * F_R * (Z0_c/(Z0_c + Z_Gn))
    
    return val

def Z_trans_along_shorted_tl(Z0, phase_vel, L, z, omega):

    tau_r = L/phase_vel
    tau_z = z/phase_vel

    val = 1j * Z0 * (np.tan(tau_r*omega)*np.cos(tau_z * omega) - np.sin(tau_z * omega)) 

    return val

def Z_trans_along_open_tl(Z0, phase_vel, L, z, omega):

    tau_r = L/phase_vel
    tau_z = z/phase_vel

    val = -1j * Z0 * (np.cos(tau_z * omega)/np.tan(tau_r*omega) + np.sin(tau_z * omega)) 

    return val

def voltage_at_source_location(Z0, phase_vel, Cm_per_len, l_c, l_Gf, l_Gn, omega):

    ### approximate solution that ignores the backaction from the coupled circuit

    Cl, Ll = transmission_line_C_and_L(phase_vel, Z0)

    Z0_c = Zchar(Cl + Cm_per_len, Ll)
    phase_vel_c = omega_r(Cl + Cm_per_len, Ll)

    tau_f = l_Gf / phase_vel

    Zinput_f = Zinput(Z0, 0, omega * tau_f)

    tau_c = l_c/phase_vel_c

    Zinput_c = Zinput(Z0_c, Zinput_f, omega * tau_c)

    tau_n = l_Gn / phase_vel

    Zinput_n = Zinput(Z0, Zinput_c, omega * tau_n)

    I_in = 1

    v_in = Zinput_n * I_in

    v_out, I_out = transmission_line_voltage_current_out(Z0, tau_n*omega, v_in, I_in)

    val = v_out

    return val

def find_notch_filter_frequency_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=3*10**8/2.5, Z0=65, search_span = 4 * 2*np.pi * 1e9, search_spacing=(10*2*np.pi*10**6), receiver_type = 'lambda/4'):
    
    # find the filter frequency by solving the analytic equation using the weak coupling assumption

    # Define the two sides of the equation
    # We will find the zero crossing, to get the solution to the equation
    # Note: omega is the variable in which we want to find the crossing

    def defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omega):
        Cl, Ll= transmission_line_C_and_L(phase_vel, Z0)
        phase_vel_c = omega_r(Cl + Cm, Ll)
        tau_c = l_c / phase_vel_c
        tau_G_dash = (l_Gf / phase_vel + l_c / phase_vel_c / 2) 
        tau_R_dash = (l_Rf / phase_vel + l_c / phase_vel_c / 2) 
        val = omega * tau_c / (np.sin(omega * tau_c)) * (np.cos(2*omega*tau_G_dash) + np.cos(2*omega*tau_R_dash)) /(1 + np.cos(2*omega*(tau_G_dash + tau_R_dash)))
        return val

    def defining_eq2(Z0, phase_vel, Lm, Cm):
        Cl, Ll= transmission_line_C_and_L(phase_vel, Z0)
        Zm = Zchar(Cm, Lm)
        Z0_c = Zchar(Cl + Cm, Ll)
        val = ((Zm/Z0_c)**2 + 1) / (1 - (Zm/Z0_c)**2)
        return val

    # find get rule of thumb notch freq first

    omega_approx = notch_filter_frequency_rule_of_thumb(l_c, l_Gf,l_Rf, Cm, phase_vel=phase_vel, Z0=Z0)

    if  receiver_type == 'lambda/2':
        omega_approx = omega_approx*2
        return omega_approx

    min_search = omega_approx - search_span/2
    max_search = omega_approx + search_span/2

    # Define vector containing results of equation for different omega values
    omegas = np.arange(min_search, max_search, search_spacing)
    results = defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas) - defining_eq2(Z0, phase_vel, Lm, Cm)

    # The equation will have a zero crossing once the sign of the difference changes
    asign = np.sign(results)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    signchange[0] = 0
    idxs = np.array(np.nonzero(signchange))[0]
    # print("idxs -> ", idxs)

    # plt.plot(omegas,  defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas) - defining_eq2(Z0, phase_vel, Lm, Cm), marker=".")
    # plt.plot(omegas[idxs], defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas[idxs]), color="red", marker="x")
    # plt.show()
    #quit()
    
    # added debugger - check for continuity of eq1 around idx:
    gaps = defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas[idxs + 1]) - defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas[idxs - 1])

    # print('gaps:', gaps)
    idx = idxs[(abs(gaps) < 30)]
    if idx.size == 0:
        print('idxs:', idxs)
        raise ValueError('No valid solution to notch frequency equation for given input parameters in specified frequency range. Therefore cannot proceed to finding Lg and Cg.')

    Z_transfer_vals = np.abs(Z_transfer_weak_coupling(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas[idx]))

    min_idx = np.argmin(Z_transfer_vals)

    idx = idx[min_idx]

    omega_f_rough  = omegas[idx]
    fine_search_spacing = 100 * 2*np.pi
    omegas = np.arange(omega_f_rough - search_spacing, omega_f_rough + search_spacing, fine_search_spacing)
    results = defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas) - defining_eq2(Z0, phase_vel, Lm, Cm)
    
    asign = np.sign(results)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    signchange[0] = 0

    idx= np.nonzero(signchange)

    val  = omegas[idx][0]

    return val

def notch_filter_frequency_rule_of_thumb(l_c, l_Gf, l_Rf, Cm = None, phase_vel=3*10**8/2.5, Z0=65):

    Cl, Ll= transmission_line_C_and_L(phase_vel, Z0)

    if Cm is not None:
        # add the effect of the scaled v_c. For weak coupling it is a small effect.

        phase_vel_c = omega_r(Cl + Cm, Ll)

    else:
        phase_vel_c = phase_vel

    tau_c = l_c / phase_vel_c
    tau_G = l_Gf / phase_vel 
    tau_R = l_Rf / phase_vel

    omega = np.pi/(2*(tau_G + tau_R + tau_c))

    return omega

def Z_input_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50, C_env = 0):

    C_val, L_val = transmission_line_C_and_L(phase_vel, Z0)

    Cm = Cm_per_len

    Lm = Lm_per_len

    L_test = np.array([[L_val, Lm],[Lm, L_val]])

    C_mutual_test = np.array([[C_val, Cm], [Cm, C_val]])

    C_Maxwell_test = C_mutual_to_Maxwell(C_mutual_test)

    #print('omegas:', omegas)
    #print('np.size(omegas):', np.size(omegas))

    if np.size(omegas) == 1 and type(omegas) is not np.array:
        omegas = np.array([omegas])

    Z_input_total_exact = []

    for omega in omegas:

        Z_env = Zline + Zcap(C_env, omega)

        Z_nb = np.diag(np.array([0, Z_input_tl(Z0, Z_env, phase_vel, l_Rn, omega)]))
        Z_fb = np.diag(np.array([Z_short_tl(Z0, phase_vel, l_Gf, omega), Z_short_tl(Z0, phase_vel, l_Rf, omega)]))

        Vs_nb = np.array([1,0])
        Vs_fb = np.array([0,0])

        V_sols_n, V_sols_f, I_sols_n, I_sols_f = V_I_solutions_sym_three(C_Maxwell_test, L_test, Z_nb, Z_fb, Vs_nb, Vs_fb, l_c, omega)

        Z_input_coupled_system = V_sols_n[0] / I_sols_n[0]

        Z_input_total_exact_val = Z_input_tl(Z0, Z_input_coupled_system, phase_vel, l_Gn, omega)

        Z_input_total_exact.append(Z_input_total_exact_val)

    Z_input_total_exact = np.array(Z_input_total_exact)

    if np.size(omegas) == 1 and type(omegas) is not np.array:
        Z_input_total_exact = Z_input_total_exact[0]

    #print('Z_input_total_exact:', Z_input_total_exact)

    return Z_input_total_exact

############################
### Lumped Element Model ###
############################

def transmission_line_C_and_L(phase_vel, Z0):

    C_val = 1/(phase_vel*Z0)
    L_val = Z0 / phase_vel

    return C_val, L_val

def transmission_line_Zchar_phase_vel(Lmutual, Cmutual):

    Zchar = np.sqrt(Lmutual/ Cmutual)

    phasevel = 1/np.sqrt(Lmutual * Cmutual)

    return Zchar, phasevel

def lumped_model_C_and_L(phase_vel, Z0, cpw__length, res_type = 'lambda/4'):

    tl_C_val, tl_L_val = transmission_line_C_and_L(phase_vel, Z0)

    if res_type == 'lambda/4':
        C_val = tl_C_val * cpw__length / 2
        L_val = 8 * tl_L_val * cpw__length / np.pi**2 

    elif res_type == 'lambda/2':
        L_val = tl_L_val * cpw__length / 2
        C_val = 2*tl_C_val * cpw__length / np.pi**2 

    return C_val, L_val

def lumped_model_Cg_and_Lg(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, receiver_type = 'lambda/4'):

    omega_f = find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0, search_span=2*2*np.pi*10**9, search_spacing=1e6*2*np.pi, receiver_type = receiver_type)
    #print('omega_f/2pi (GHz):', omega_f/(2*np.pi*1e9))
    Z0_f = find_notch_filter_char_impedance(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega_f, phase_vel=phase_vel, Z0=Z0, receiver_type = receiver_type)
    #print('Z0_f:', Z0_f)

    Cg = 1/(omega_f*Z0_f)
    Lg =  Z0_f / omega_f

    return Cg, Lg

def get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, receiver_type = 'lambda/4'):
        
    cpw__length1 = l_c + l_Gf + l_Gn
    cpw__length2 = l_c + l_Rf + l_Rn
    
    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1)
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2, res_type = receiver_type)
    Cg, Lg = lumped_model_Cg_and_Lg(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, receiver_type = receiver_type)
    
    # Cg_basic = Cm_per_len *l_c
    # Lg_basic = L1*L2/ (Lm_per_len * l_c)

    return C1, L1, C2, L2, Cg, Lg

def lumped_model_Z_transmission(omega, C1, L1, C2, L2, Cg, Lg, Cg2 = None, Lg2 = None):

    Z_res1 = Zres(C1, L1, omega)
    Z_res_coupler = Zres(Cg, Lg, omega)
    if Cg2 is not None and Lg2 is not None:
        Z_res_coupler = Z_res_coupler + Zres(Cg2, Lg2, omega)
    Z_res2 = Zres(C2, L2, omega)

    val = Z_res1 * Z_res2 / Z_res_coupler * 1/(1 + (Z_res1 + Z_res2)/Z_res_coupler)

    return val 

def lumped_model_resonator_coupling(C1, L1, C2, L2, Cg, Lg):

    omega_1 = omega_r(C1, L1)
    omega_2 = omega_r(C2, L2)

    Z21_omega_1 = lumped_model_Z_transmission(omega_1, C1, 1e3, C2, 1e3, Cg, Lg)
    Z21_omega_1 = lumped_model_Z_transmission(omega_2, C1, 1e3, C2, 1e3, Cg, Lg)

    val = -0.25*(omega_1**3*omega_2**3*C1*C2)**0.5*(np.imag(Z21_omega_1)/omega_1 + np.imag(Z21_omega_1)/omega_2)

    return val

def find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=3*10**8/2.5, Z0=65, search_span=2*2*np.pi*10**9, search_spacing=(2.5*2*np.pi*10**6), receiver_type = 'lambda/4'):

    ## find notch frequency by directly solving for the zeros of the transfer impedance

    # Define vector containing results of equation for different omega values

    omega_approx = find_notch_filter_frequency_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=phase_vel, Z0=Z0, search_span = search_span * 2, search_spacing=search_spacing * 2, receiver_type = receiver_type)

    min_search = omega_approx - search_span/2
    max_search = omega_approx + search_span/2

    omegas = np.arange(min_search, max_search, search_spacing)

    Z_vals = np.imag(Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0, receiver_type = receiver_type))

    # The transmission will have a zero crossing at the notch
    asign = np.sign(Z_vals)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    signchange[0] = 0
    idxs = np.array(np.nonzero(signchange))[0]

    # plt.plot(omegas/(2*np.pi*1e9), Z_vals)
    # plt.plot(omegas[idxs]/(2*np.pi*1e9),[0]*len(idxs), marker = 'x')
    # plt.show()

    # print("idxs -> ", idxs)
    
    # added debugger - check for continuity of transmission around idx:
    gaps = Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas[idxs + 1], phase_vel=phase_vel, Z0=Z0, receiver_type = receiver_type) -  Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas[idxs - 1], phase_vel=phase_vel, Z0=Z0, receiver_type = receiver_type)
    
    gaps = gaps.ravel()

    #print('idxs:', idxs)
    #print('gaps:', gaps)
    
    idx = idxs[(abs(gaps) < 10000)]
    if idx.size == 0:
        print('idxs:', idxs)
        raise ValueError('No valid solution to notch frequency equation for given input parameters in specified frequency range. Therefore cannot proceed to finding Lg and Cg.')

    val  = omegas[idx][0]

    return val

def Z_transfer_differential(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f, phase_vel=3*10**8/2.5, Z0=65, delta_omega = 15 * 2*np.pi, receiver_type = 'lambda/4'):

    # delta_omega is 15 Hz by default

    Z_transfer_plus = Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f + delta_omega/2, phase_vel=phase_vel, Z0=Z0, receiver_type = receiver_type)
    Z_transfer_minus = Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f - delta_omega/2, phase_vel=phase_vel, Z0=Z0, receiver_type = receiver_type)

    grad_Ztransfer = (Z_transfer_plus - Z_transfer_minus)/(delta_omega)

    grad_Ztransfer = grad_Ztransfer[0]

    return grad_Ztransfer

def find_notch_filter_char_impedance(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f, phase_vel=3*10**8/2.5, Z0=65, receiver_type = 'lambda/4'):

    cpw__length1 = l_c + l_Gf + l_Gn
    cpw__length2 = l_c + l_Rf + l_Rn

    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1)
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2, res_type=receiver_type)

    Zres_1_omega_f = Zres(C1, L1, omega_f)
    Zres_2_omega_f = Zres(C2, L2, omega_f)

    grad_Ztransfer_val = Z_transfer_differential(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f, phase_vel=phase_vel, Z0=Z0, receiver_type = receiver_type)

    #print('grad_Ztransfer_val:', grad_Ztransfer_val)

    val = np.real(1j*2*Zres_1_omega_f * Zres_2_omega_f / (omega_f * grad_Ztransfer_val))

    return val

def Z_notch_symbolic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65):
    
    tl_C, tl_L_val = transmission_line_C_and_L(phase_vel, Z0)

    cpw__length1 = l_c + l_Gf + l_Gn
    cpw__length2 = l_c + l_Rf + l_Rn
    
    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1)
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2)

    omega_n = notch_filter_frequency_rule_of_thumb(l_c, l_Gf, l_Rf, phase_vel=phase_vel, Z0=Z0)

    omega_r = 1/np.sqrt(C1 * L1)
    omega_p = 1/np.sqrt(C2 * L2)

    val_1 = Z0 * (64/np.pi**3) * np.cos(omega_n * np.pi/(2*omega_r)) * np.cos(omega_n * np.pi/(2*omega_p)) / (((omega_r/omega_n)-(omega_n/omega_r))*((omega_p/omega_n)-(omega_n/omega_p)))

    tau_c = l_c/phase_vel

    val_2 = (tl_C/Cm)/np.sin(omega_n * tau_c)

    val = val_1*val_2

    return val

def get_lumped_elements_from_symbolic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, receiver_type = 'lambda/4'):
        
    cpw__length1 = l_c + l_Gf + l_Gn
    cpw__length2 = l_c + l_Rf + l_Rn
    
    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1)
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2, res_type = receiver_type)
    #Cg, Lg = lumped_model_Cg_and_Lg(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, receiver_type = receiver_type)
    
    omega_notch_symbolic = notch_filter_frequency_rule_of_thumb(l_c, l_Gf, l_Rf, phase_vel=phase_vel, Z0=Z0)

    Z_notch_symbolic_val = Z_notch_symbolic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Cm_per_len, phase_vel=phase_vel, Z0=Z0)

    Cg = 1/(Z_notch_symbolic_val*omega_notch_symbolic)
    Lg = Z_notch_symbolic_val/omega_notch_symbolic

    # Cg_basic = Cm_per_len *l_c
    # Lg_basic = L1*L2/ (Lm_per_len * l_c)

    return C1, L1, C2, L2, Cg, Lg

# lines coupled by LE capacitor & inductor

def get_lumped_elements_LE_coupling(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, phase_vel=3*10**8/2.5, Z0=65):

    cpw__length1 = l_Gf + l_Gn
    cpw__length2 = l_Rf + l_Rn
    
    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1)
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2)
    Cg, Lg = lumped_model_Cg_and_Lg_LE_coupling(phase_vel, Z0, l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc)

    return C1, L1, C2, L2, Cg, Lg

def lumped_model_Cg_and_Lg_LE_coupling(phase_vel, Z0, l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc):

    omega_f = find_notch_filter_frequency_LE_coupling(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, phase_vel, Z0, 1e9*2*np.pi, 10e9*2*np.pi, 1e6*2*np.pi)
    #print('omega_f/2pi (GHz):', omega_f/(2*np.pi*1e9))
    Z0_f = find_notch_filter_char_impedance_LE_coupling(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, omega_f, phase_vel=phase_vel, Z0=Z0)
    #print('Z0_f:', Z0_f)

    Cg = 1/(omega_f*Z0_f)
    Lg =  Z0_f / omega_f

    return Cg, Lg

def find_notch_filter_frequency_LE_coupling(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, phase_vel=3*10**8/2.5, Z0=65, min_search=3*2*np.pi*10**9, max_search=13*2*np.pi*10**9, search_spacing=(10*2*np.pi*10**6)):

    ## find notch frequency by directly solving for the zeros of the transfer impedance

    # Define vector containing results of equation for different omega values
    omegas = np.arange(min_search, max_search, search_spacing)

    Z_vals = np.imag(Z_transfer_LE_coupled_lines(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, omegas, phase_vel=phase_vel, Z0=Z0))

    # The transmission will have a zero crossing at the notch
    asign = np.sign(Z_vals)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    signchange[0] = 0
    idxs = np.array(np.nonzero(signchange))[0]

    # plt.plot(omegas/(2*np.pi*1e9), Z_vals)
    # plt.plot(omegas[idxs]/(2*np.pi*1e9),[0]*len(idxs), marker = 'x')
    # plt.show()

    # print("idxs -> ", idxs)

    #quit()
    
    # added debugger - check for continuity of transmission around idx:
    gaps = Z_transfer_LE_coupled_lines(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, omegas[idxs + 1], phase_vel=phase_vel, Z0=Z0) -  Z_transfer_LE_coupled_lines(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, omegas[idxs - 1], phase_vel=phase_vel, Z0=Z0)
    
    #print('idxs:', idxs)
    #print('gaps:', gaps)
    
    idx = idxs[(abs(gaps) < 30)]
    if idx.size == 0:
        print('idxs:', idxs)
        raise ValueError('No valid solution to notch frequency equation for given input parameters in specified frequency range. Therefore cannot proceed to finding Lg and Cg.')

    val  = omegas[idx][0]

    return val

def find_notch_filter_char_impedance_LE_coupling(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, omega_f, phase_vel=3*10**8/2.5, Z0=65):

    cpw__length1 =  l_Gf + l_Gn
    cpw__length2 =  l_Rf + l_Rn

    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1)
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2)

    Zres_1_omega_f = Zres(C1, L1, omega_f)
    Zres_2_omega_f = Zres(C2, L2, omega_f)

    delta_omega = 10 * 2*np.pi # 10 Hz

    Z_transfer_plus = Z_transfer_LE_coupled_lines(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, omega_f + delta_omega/2, phase_vel=phase_vel, Z0=Z0)
    Z_transfer_minus = Z_transfer_LE_coupled_lines(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, omega_f - delta_omega/2, phase_vel=phase_vel, Z0=Z0)

    grad_Ztransfer = (Z_transfer_plus - Z_transfer_minus)/(delta_omega)

    val = np.real(1j*2*Zres_1_omega_f * Zres_2_omega_f / (omega_f * grad_Ztransfer))

    return val

###################################
### Z matrix transfer functions ###
###################################

def Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=3*10**8/2.5, Z0=65, receiver_type = 'lambda/4', receiver_side = 'near'):

    C_val, L_val = transmission_line_C_and_L(phase_vel, Z0)

    Cm = Cm_per_len

    Lm = Lm_per_len

    C_mutual_test = np.array([[C_val, Cm], [Cm, C_val]])

    C_Maxwell_test = C_mutual_to_Maxwell(C_mutual_test)

    L_test = np.array([[L_val, Lm],[Lm, L_val]])

    V_sols_n_arr = []
    V_sols_f_arr = []
    I_sols_n_arr = []
    I_sols_f_arr = []
    V_from_input_I_arr = []
    tl_voltage_scale_factor_arr = []
    Z_input_exact_arr = []

    # if type(omegas) is not list() or type(omegas) is not np.array():
    #     omegas = [omegas]

    #print(np.size(omegas))
    if np.size(omegas) ==1:
        omegas = [omegas]

    # if len(omegas) == 0:
    #     #if type(omegas) != list() or type(omegas) != np.array():
    #     omegas = [omegas]

    for omega in omegas:

        Z_nb = np.diag([0, complex(Z_open_tl(Z0, phase_vel, l_Rn, omega))]) ## np.array(

        #Z_nb = np.asarray(Z_nb,dtype=object)

        if receiver_type == 'lambda/4':
            Z_fb = np.diag(np.array([Z_short_tl(Z0, phase_vel, l_Gf, omega), Z_short_tl(Z0, phase_vel, l_Rf, omega)]))
        if receiver_type == 'lambda/2':
            Z_fb = np.diag(np.array([Z_short_tl(Z0, phase_vel, l_Gf, omega), Z_open_tl(Z0, phase_vel, l_Rf, omega)]))

        Vs_nb = np.array([1,0])
        Vs_fb = np.array([0,0])

        V_sols_n, V_sols_f, I_sols_n, I_sols_f = V_I_solutions_sym_three(C_Maxwell_test, L_test, Z_nb, Z_fb, Vs_nb, Vs_fb, l_c, omega)

        V_sols_n_arr.append(V_sols_n)
        V_sols_f_arr.append(V_sols_f)

        I_sols_n_arr.append(I_sols_n)
        I_sols_f_arr.append(I_sols_f)

        Z_input_exact = V_sols_n[0] / I_sols_n[0]

        Z_input_exact_arr.append(Z_input_exact)

        V_from_input_I = voltage_at_source_location_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, Z0, phase_vel, omega, solved_Z_input_coupled_system=Z_input_exact)

        V_from_input_I_arr.append(V_from_input_I)

        tau_Rn = l_Rn/phase_vel
        tau_Rf = l_Rf/phase_vel

        if receiver_side == 'near':

            tl_voltage_scale_factor = 1/np.cos(tau_Rn*omega)

        if receiver_side == 'far':
        
            tl_voltage_scale_factor = 1/np.cos(tau_Rf*omega)

        tl_voltage_scale_factor_arr.append(tl_voltage_scale_factor)

    V_sols_n_arr = np.array(V_sols_n_arr)    
    V_sols_f_arr = np.array(V_sols_f_arr)
    V_from_input_I_arr = np.array(V_from_input_I_arr)
    tl_voltage_scale_factor_arr = np.array(tl_voltage_scale_factor_arr)

    Z_input_exact_arr = np.array(Z_input_exact_arr)

    ## take the imag part to remove tiny numerical errors

    if receiver_side == 'near':
        Z_transfer_exact = 1j*np.imag(tl_voltage_scale_factor_arr * V_from_input_I_arr * V_sols_n_arr[:, -1])

    if receiver_side == 'far':
        Z_transfer_exact = 1j*np.imag(tl_voltage_scale_factor_arr * V_from_input_I_arr * V_sols_f_arr[:, -1])

    return Z_transfer_exact

def Z_transfer_weak_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega, phase_vel=3*10**8/2.5, Z0=65):

    L_G = l_Gn + l_c + l_Gf

    I_dummy_in = 1

    #L_G = l_c + l_Gf

    ## approximate solution
    #Z_trans_along_shorted_tr_val = Z_trans_along_shorted_tl(Z0, phase_vel, L_G, l_Gn, omega)

    ## more accurate solution
    Z_trans_along_shorted_tr_val = voltage_at_source_location(Z0, phase_vel, Cm_per_len, l_c, l_Gf, l_Gn, omega) / I_dummy_in

    voltage_transmission_coupled_lines_val = voltage_transmission_coupled_lines(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega)

    tau_Rn = l_Rn/phase_vel

    Rn_tl_voltage_scale_factor = 1/np.cos(tau_Rn*omega)

    val = Z_trans_along_shorted_tr_val * voltage_transmission_coupled_lines_val * Rn_tl_voltage_scale_factor

    return val

def Z_transfer_equivalent_LE_circuit(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega, phase_vel=3*10**8/2.5, Z0=65, receiver_type = 'lambda/4'):

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0, receiver_type = receiver_type)

    val = lumped_model_Z_transmission(omega, C1, L1, C2, L2, Cg, Lg)

    return val

def Z_transfer_equivalent_LE_circuit_from_symbolic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Cm_per_len, omega, phase_vel=3*10**8/2.5, Z0=65, receiver_type = 'lambda/4'):

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements_from_symbolic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Cm_per_len, phase_vel, Z0, receiver_type = receiver_type)

    val = lumped_model_Z_transmission(omega, C1, L1, C2, L2, Cg, Lg)

    return val

def Z_transfer_LE_coupled_lines(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, omegas, phase_vel=3*10**8/2.5, Z0=65):

    Z_Gf = Z_short_tl(Z0, phase_vel, l_Gf, omegas)
    Z_Rf = Z_short_tl(Z0, phase_vel, l_Rf, omegas)
    Z_Rn = Z_open_tl(Z0, phase_vel, l_Rn, omegas)

    Z_LE_coupler = Zpara(Zind(Lc, omegas), Zcap(Cc, omegas))

    Z_load = Zpara(Z_Gf, Z_LE_coupler + Zpara(Z_Rf, Z_Rn))

    Z_in = Z_input_tl(Z0, Z_load, phase_vel, l_Gn, omegas)

    I_in = 1

    V_in = Z_in * I_in

    electric_length_Gn = l_Gn/phase_vel * omegas

    V_out, I_out = transmission_line_voltage_current_out(Z0, electric_length_Gn, V_in, I_in)

    I_down = I_out * Z_load /  (Z_LE_coupler + Zpara(Z_Rf, Z_Rn))

    V_in2 = V_out - I_down * Z_LE_coupler

    I_in2 = I_down * Zpara(Z_Rf, Z_Rn) / Z_Rn

    electric_length_Rn = l_Rn/phase_vel * omegas

    V_out2 = V_in2/np.cos(electric_length_Rn)

    val = V_out2 / I_in

    return val

def Z_transfer_equivalent_LE_circuit_LE_coupling(l_Gf, l_Gn, l_Rf, l_Rn,  Lc, Cc, omega, phase_vel=3*10**8/2.5, Z0=65):

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements_LE_coupling(l_Gf, l_Gn, l_Rf, l_Rn, Lc, Cc, phase_vel=phase_vel, Z0=Z0)

    val = lumped_model_Z_transmission(omega, C1, L1, C2, L2, Cg, Lg)

    return val

## S matrix

def S_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=3*10**8/2.5, Z0=65, Zport = 50, receiver_type = 'lambda/4'):

    Ztrans = Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=phase_vel, Z0=Z0, receiver_type = receiver_type)

    Zinput = []

    for omega in omegas:

        Zin = Z_input_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega, phase_vel=3*10**8/2.5, Z0=65, Zline = 50, C_env = 0)

        Zinput.append(Zin)

    Zinput = np.array(Zinput)

    I_mat = np.diag([1,1])
    sqrt_Z_char_mat = np.diag([np.sqrt(Zport),np.sqrt(Zport)])
    sqrt_Y_char_mat = np.diag([1/np.sqrt(Zport),1/np.sqrt(Zport)])

    S_mat = []

    for Zin, Ztr in zip(Zinput,Ztrans):

        Z_mat = np.array([[Zin, Ztr],[Ztr, Zin]])
        S_m =  np.linalg.inv(sqrt_Y_char_mat@Z_mat@sqrt_Y_char_mat + I_mat)@(sqrt_Y_char_mat@Z_mat@sqrt_Y_char_mat - I_mat)

        S_mat.append(S_m)

    S_mat = np.array(S_mat)

    return S_mat

#################################################
###      Distributed Element Calculations     ###
### 2 capacitively coupled transmission lines ###
################ both lamdba/4 ##################
#################################################

### For 2 lamdbda/4 transmission lines coupled by a lumped element capacitor

def Z_transfer_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega, phase_vel=3*10**8/2.5, Z0=65):

    I_dummy = 1

    Z_Gf = Z_short_tl(Z0, phase_vel, l_Gf, omega)

    Z_cap = Zcap(Cm, omega)

    Z_Rn = Z_open_tl(Z0, phase_vel, l_Rn, omega)
    Z_Rf = Z_short_tl(Z0, phase_vel, l_Rf, omega) #Z_short_tl(Z0, phase_vel, l_Rf, omega)

    Z_2 = Z_cap + Zpara(Z_Rn, Z_Rf)

    Ztot = Zpara(Z_Gf, Z_2)

    ## approximate solution
    #Z_trans_along_shorted_tr_val = Z_trans_along_shorted_tl(Z0, phase_vel, L_G, l_Gn, omega)

    ## more accurate solution
    Ztot_in = Zinput(Z0, Ztot, omega*l_Gn/phase_vel)

    V_in = Ztot_in * I_dummy

    V_out, I_out = transmission_line_voltage_current_out(Z0, omega*l_Gn/phase_vel, V_in, I_dummy)

    #V0 = I_dummy * Ztot

    I2 = I_out * Ztot/Z_2

    V1 = V_out - Z_cap*I2

    #I_3 = I2*Zpara(Z_Rn, Z_Rf)/Z_Rn

    tau_Rn = l_Rn/phase_vel

    Rn_tl_voltage_scale_factor = 1/np.cos(tau_Rn*omega)

    Vout = Rn_tl_voltage_scale_factor * V1

    return Vout/I_dummy

def Z_transfer_differential_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f, phase_vel=3*10**8/2.5, Z0=65, delta_omega = 10 * 2*np.pi):

    # delta_omega is 15 Hz by default

    Z_transfer_0 = Z_transfer_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f+delta_omega/200, phase_vel=phase_vel, Z0=Z0)
    Z_transfer_1 = Z_transfer_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f + delta_omega, phase_vel=phase_vel, Z0=Z0)
    Z_transfer_2 = Z_transfer_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f + delta_omega*2, phase_vel=phase_vel, Z0=Z0)
    Z_transfer_3 = Z_transfer_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f + delta_omega*3, phase_vel=phase_vel, Z0=Z0)

    grad_Ztransfer_0 = (Z_transfer_1 - Z_transfer_0) / delta_omega
    grad_Ztransfer_1 = (Z_transfer_2 - Z_transfer_1) / delta_omega
    grad_Ztransfer_2 = (Z_transfer_3 - Z_transfer_2) / delta_omega

    grad_grad_Ztransfer_0 = (grad_Ztransfer_1 - grad_Ztransfer_0) / delta_omega
    grad_grad_Ztransfer_1 = (grad_Ztransfer_2 - grad_Ztransfer_1) / delta_omega

    grad_grad_grad_Ztransfer = (grad_grad_Ztransfer_1- grad_grad_Ztransfer_0) / delta_omega

    return grad_grad_grad_Ztransfer

def find_coupling_cap_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65):

    cpw__length1 = l_Gf + l_Gn
    cpw__length2 = l_Rf + l_Rn

    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1, res_type = 'lambda/4')
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2, res_type = 'lambda/4')

    omega_1 = np.sqrt(1/(L1*C1))
    omega_2 = np.sqrt(1/(L2*C2))

    omega_bar = (omega_1 + omega_2) / 2

    # Zchar_1 = np.sqrt(L1/C1)
    # Zchar_2 = np.sqrt(L2/C2)

    #grad_grad_grad_Ztransfer_val = Z_transfer_differential_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, Cm, 0, phase_vel=phase_vel, Z0=Z0)

    #val = 1j * omega_1 * omega_2 / (Zchar_1 * Zchar_2) * grad_grad_grad_Ztransfer_val / 6

    omega_test = 5000 * 2*np.pi *1e6

    Z21 = Z_transfer_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_test, phase_vel=phase_vel, Z0=Z0)

    ZR1 = Zres(C1, L1, omega_test)
    ZR2 = Zres(C2, L2, omega_test)

    val = -1j * Z21 / ( omega_test * ZR1*ZR2 )

    ## testing simple version heres

    omega_n_r = np.pi*phase_vel/l_Gf
    omega_n_p = np.pi*phase_vel/l_Rf

    Cg = Cm * np.sin(np.pi*omega_1/omega_n_r)* np.sin(np.pi*omega_2/omega_n_p)
    
    #Cg = Cm * np.sin(np.pi*omega_bar/omega_n_r)* np.sin(np.pi*omega_bar/omega_n_p)

    val = Cg

    return val

def lumped_model_Cg_and_Lg_direct_cap(phase_vel, Z0, l_Gf, l_Gn, l_Rf, l_Rn, Cm):

    # omega_f1 = np.pi*phase_vel / l_Gf
    # omega_f2 = np.pi*phase_vel / l_Rf

    Cc = find_coupling_cap_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=phase_vel, Z0=Z0)

    Ln = 1e3 

    return Cc, Ln

def get_lumped_elements_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65):
        
    cpw__length1 = l_Gf + l_Gn
    cpw__length2 = l_Rf + l_Rn
    
    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1, res_type = 'lambda/4')
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2, res_type = 'lambda/4')
    Cn, Ln = lumped_model_Cg_and_Lg_direct_cap(phase_vel, Z0, l_Gf, l_Gn, l_Rf, l_Rn, Cm)

    return C1, L1, C2, L2, Cn, Ln

def J_coupling_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65):

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel, Z0)
    J = lumped_model_resonator_coupling(C1, L1, C2, L2, Cg, Lg)

    return J

def J_coupling_direct_cap_analytic(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65):

    cpw__length1 = l_Gf + l_Gn
    cpw__length2 = l_Rf + l_Rn
    
    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1, res_type = 'lambda/4')
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2, res_type = 'lambda/4')

    omega_1 = 1/np.sqrt(C1 * L1)
    omega_2 = 1/np.sqrt(C2 * L2)

    omega_bar = (omega_1 + omega_2)/2

    omega_n_r = np.pi*phase_vel/l_Gf
    omega_n_p = np.pi*phase_vel/l_Rf

    Cg = Cm * np.sin(np.pi*omega_1/omega_n_r)* np.sin(np.pi*omega_2/omega_n_p)

    val = 2/np.pi * Z0 * omega_bar**2 * Cg

    return val

#################################################
###      Distributed Element Calculations     ###
### 2 capacitively coupled transmission lines ###
########## one lamdba/4 one lambda/2 ############
#################################################

# ### For one lamdbda/2 and one lambda/4 transmission lines coupled by a lumped element capacitor

def Z_transfer_direct_cap_exact_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega, phase_vel=3*10**8/2.5, Z0=65):

    I_dummy = 1

    Z_Gf = Z_short_tl(Z0, phase_vel, l_Gf, omega)

    Z_cap = Zcap(Cm, omega)

    Z_Rn = Z_open_tl(Z0, phase_vel, l_Rn, omega)
    Z_Rf = Z_open_tl(Z0, phase_vel, l_Rf, omega) #Z_short_tl(Z0, phase_vel, l_Rf, omega)

    Z_2 = Z_cap + Zpara(Z_Rn, Z_Rf)

    Ztot = Zpara(Z_Gf, Z_2)

    ## approximate solution
    #Z_trans_along_shorted_tr_val = Z_trans_along_shorted_tl(Z0, phase_vel, L_G, l_Gn, omega)

    ## more accurate solution
    Ztot_in = Zinput(Z0, Ztot, omega*l_Gn/phase_vel)

    V_in = Ztot_in * I_dummy

    V_out, I_out = transmission_line_voltage_current_out(Z0, omega*l_Gn/phase_vel, V_in, I_dummy)

    #V0 = I_dummy * Ztot

    I2 = I_out * Ztot/Z_2

    V1 = V_out - Z_cap*I2

    #I_3 = I2*Zpara(Z_Rn, Z_Rf)/Z_Rn

    tau_Rn = l_Rn/phase_vel

    Rn_tl_voltage_scale_factor = 1/np.cos(tau_Rn*omega)

    Vout = Rn_tl_voltage_scale_factor * V1

    return Vout/I_dummy

def Z_transfer_differential_direct_cap_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f, phase_vel=3*10**8/2.5, Z0=65, delta_omega = 15 * 2*np.pi):

    # delta_omega is 15 Hz by default

    Z_transfer_0 = Z_transfer_direct_cap_exact_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f-delta_omega/2, phase_vel=phase_vel, Z0=Z0)
    Z_transfer_1 = Z_transfer_direct_cap_exact_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f+delta_omega/2, phase_vel=phase_vel, Z0=Z0)

    grad_Ztransfer = (Z_transfer_1 - Z_transfer_0) / delta_omega

    return grad_Ztransfer

def find_notch_filter_char_impedance_direct_cap_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f, phase_vel=3*10**8/2.5, Z0=65):

    cpw__length1 = l_Gf + l_Gn
    cpw__length2 = l_Rf + l_Rn

    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1, res_type='lambda/4')
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2, res_type='lambda/2')

    Zres_1_omega_f = Zres(C1, L1, omega_f)
    Zres_2_omega_f = Zres(C2, L2, omega_f)

    grad_Ztransfer_val = Z_transfer_differential_direct_cap_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f, phase_vel=phase_vel, Z0=Z0)

    #print('grad_Ztransfer_val:', grad_Ztransfer_val)

    val = np.real(1j*2*Zres_1_omega_f * Zres_2_omega_f / (omega_f * grad_Ztransfer_val))

    return val

def lumped_model_Cg_and_Lg_direct_cap_B(phase_vel, Z0, l_Gf, l_Gn, l_Rf, l_Rn, Cm):

    omega_f = np.pi*phase_vel / (2*l_Rf)
    #omega_f2 = np.pi*phase_vel / l_Rf

    #print('omega_f(GHz):', omega_f / (2*np.pi*1e9))

    #print('omega_f/2pi (GHz):', omega_f/(2*np.pi*1e9))
    Z0_f = find_notch_filter_char_impedance_direct_cap_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f, phase_vel=phase_vel, Z0=Z0)
    #print('Z0_f:', Z0_f)

    Cg = 1/(omega_f*Z0_f)
    Lg =  Z0_f / omega_f

    return Cg, Lg

def get_lumped_elements_direct_cap_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65):
        
    cpw__length1 = l_Gf + l_Gn
    cpw__length2 = l_Rf + l_Rn
    
    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1, res_type = 'lambda/4')
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2, res_type = 'lambda/2')
    Cn, Ln = lumped_model_Cg_and_Lg_direct_cap_B(phase_vel, Z0, l_Gf, l_Gn, l_Rf, l_Rn, Cm)

    return C1, L1, C2, L2, Cn, Ln

def J_coupling_direct_cap_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65):

    l_Rf_for_J = max(l_Rf, l_Rn)
    l_Rn_for_J = min(l_Rf, l_Rn)

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements_direct_cap_B(l_Gf, l_Gn, l_Rf_for_J, l_Rn_for_J, Cm, phase_vel, Z0)
    J = lumped_model_resonator_coupling(C1, L1, C2, L2, Cg, Lg)

    return J

def J_coupling_direct_cap_analytic_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65):

    cpw__length1 = l_Gf + l_Gn
    cpw__length2 = l_Rf + l_Rn
    
    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1, res_type = 'lambda/4')
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2, res_type = 'lambda/2')

    omega_r = 1/np.sqrt(C1 * L1)
    omega_p = 1/np.sqrt(C2 * L2)

    omega_bar = (omega_r + omega_p) / 2

    tau_Rf = l_Rf/phase_vel
    tau_Rn = l_Rn/phase_vel

    tau_Gf = l_Gf/phase_vel

    tau_notch = max(tau_Rf, tau_Rn)

    omega_n = np.pi / (2 * tau_notch)

    #print('omega_n:', omega_n / (2*np.pi*1e9))

    J_val = np.pi **2 / 32 * Z0 * Cm * omega_n * omega_bar * (omega_r/omega_n - omega_n/omega_r) * (omega_p/omega_n - omega_n/omega_p) * (omega_bar/omega_n - omega_n/omega_bar) * np.sin(omega_n * tau_Gf) * 1/(np.cos(omega_n * np.pi / (2*omega_r)) * np.sin(omega_n * np.pi / (omega_p)))

    return J_val

#################################################
###      Distributed Element Calculations     ###
### 2 capacitively coupled transmission lines ###
################ both lamdba/2 ##################
#################################################

# ### For two lamdbda/2 transmission lines coupled by a lumped element capacitor

def Z_transfer_direct_cap_exact_C(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega, phase_vel=3*10**8/2.5, Z0=65):

    I_dummy = 1

    Z_Gf = Z_open_tl(Z0, phase_vel, l_Gf, omega) #Z_short_tl(Z0, phase_vel, l_Gf, omega)

    Z_cap = Zcap(Cm, omega)

    Z_Rn = Z_open_tl(Z0, phase_vel, l_Rn, omega)
    Z_Rf = Z_open_tl(Z0, phase_vel, l_Rf, omega) #Z_short_tl(Z0, phase_vel, l_Rf, omega)

    Z_2 = Z_cap + Zpara(Z_Rn, Z_Rf)

    Ztot = Zpara(Z_Gf, Z_2)

    ## approximate solution
    #Z_trans_along_shorted_tr_val = Z_trans_along_shorted_tl(Z0, phase_vel, L_G, l_Gn, omega)

    ## more accurate solution
    Ztot_in = Zinput(Z0, Ztot, omega*l_Gn/phase_vel)

    V_in = Ztot_in * I_dummy

    V_out, I_out = transmission_line_voltage_current_out(Z0, omega*l_Gn/phase_vel, V_in, I_dummy)

    #V0 = I_dummy * Ztot

    I2 = I_out * Ztot/Z_2

    V1 = V_out - Z_cap*I2

    #I_3 = I2*Zpara(Z_Rn, Z_Rf)/Z_Rn

    tau_Rn = l_Rn/phase_vel

    Rn_tl_voltage_scale_factor = 1/np.cos(tau_Rn*omega)

    Vout = Rn_tl_voltage_scale_factor * V1

    return Vout/I_dummy

def Z_transfer_differential_direct_cap_C(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f, phase_vel=3*10**8/2.5, Z0=65, delta_omega = 15 * 2*np.pi):

    # delta_omega is 15 Hz by default

    Z_transfer_0 = Z_transfer_direct_cap_exact_C(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f-delta_omega/2, phase_vel=phase_vel, Z0=Z0)
    Z_transfer_1 = Z_transfer_direct_cap_exact_C(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f+delta_omega/2, phase_vel=phase_vel, Z0=Z0)

    grad_Ztransfer = (Z_transfer_1 - Z_transfer_0) / delta_omega

    return grad_Ztransfer

def find_notch_filter_char_impedance_direct_cap_C(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f, phase_vel=3*10**8/2.5, Z0=65):

    cpw__length1 = l_Gf + l_Gn
    cpw__length2 = l_Rf + l_Rn

    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1, res_type='lambda/4')
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2, res_type='lambda/2')

    Zres_1_omega_f = Zres(C1, L1, omega_f)
    Zres_2_omega_f = Zres(C2, L2, omega_f)

    grad_Ztransfer_val = Z_transfer_differential_direct_cap_C(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f, phase_vel=phase_vel, Z0=Z0)

    #print('grad_Ztransfer_val:', grad_Ztransfer_val)

    val = np.real(1j*2*Zres_1_omega_f * Zres_2_omega_f / (omega_f * grad_Ztransfer_val))

    return val

def lumped_model_Cg_and_Lg_direct_cap_C(phase_vel, Z0, l_Gf, l_Gn, l_Rf, l_Rn, Cm):

    omega_f1 = np.pi*phase_vel / (2*l_Rf)
    omega_f2 = np.pi*phase_vel / (2*l_Gf)

    #print('omega_f1(GHz):', omega_f1 / (2*np.pi*1e9))
    #print('omega_f2(GHz):', omega_f2 / (2*np.pi*1e9))

    #print('omega_f/2pi (GHz):', omega_f/(2*np.pi*1e9))
    Z0_f1 = find_notch_filter_char_impedance_direct_cap_C(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f1, phase_vel=phase_vel, Z0=Z0)
    Z0_f2 = find_notch_filter_char_impedance_direct_cap_C(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omega_f2, phase_vel=phase_vel, Z0=Z0)

    #print('Z0_f:', Z0_f)

    Cg1 = 1/(omega_f1*Z0_f1)
    Lg1 =  Z0_f1 / omega_f1

    Cg2 = 1/(omega_f2*Z0_f2)
    Lg2 =  Z0_f2 / omega_f2

    return Cg1, Lg1, Cg2, Lg2

def get_lumped_elements_direct_cap_C(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65):
        
    cpw__length1 = l_Gf + l_Gn
    cpw__length2 = l_Rf + l_Rn
    
    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1, res_type = 'lambda/2')
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2, res_type = 'lambda/2')
    Cn1, Ln1, Cn2, Ln2 = lumped_model_Cg_and_Lg_direct_cap_C(phase_vel, Z0, l_Gf, l_Gn, l_Rf, l_Rn, Cm)

    return C1, L1, C2, L2, Cn1, Ln1, Cn2, Ln2

#######################################
### qubit radiative decay functions ###
#######################################

def qubit_radiative_decay_sym_3_lines_exact(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50):

    Z_input_exact = Z_input_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=phase_vel, Z0=Z0, Zline = Zline, C_env = C_ext)

    Zq_env = Zcap(C_g, omegas) + Z_input_exact

    T1 = qubit_radiative_T1(C_q, Zq_env)

    return T1

def qubit_radiative_decay_equivalent_LE_circuit(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50):

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0)

    Z1 = 1/(1j*omegas*C1 + 1/(1j*omegas*L1))
    Z2 = 1/(1j*omegas*Cg + 1/(1j*omegas*Lg))
    Z3 = 1/(1j*omegas*C2 + 1/(1j*omegas*L2))

    Zg = Zcap(C_g, omegas)
    Z_ext = Zcap(C_ext, omegas) + Zline 

    Zq_env = Zg + Zpara(Z1, Z2 + Zpara(Z3, Z_ext)) 

    T1 = qubit_radiative_T1(C_q, Zq_env)

    return T1

def qubit_radiative_decay_equivalent_LE_circuit_without_notch(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50):
    
    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0)

    J_val = lumped_model_resonator_coupling(C1, L1, C2, L2, Cg, Lg)

    eff_Cg = get_eff_Cf_from_J(l_c, l_Gf, l_Gn, l_Rf, l_Rn, J_val, phase_vel, Z0)

    #print('eff_Cg:', eff_Cg)

    Z1 = 1/(1j*omegas*C1 + 1/(1j*omegas*L1))
    Z2 = 1/(1j*omegas*eff_Cg)
    Z3 = 1/(1j*omegas*C2 + 1/(1j*omegas*L2))

    Zg = Zcap(C_g, omegas)
    Z_ext = Zcap(C_ext, omegas) + Zline 

    Zq_env = Zg + Zpara(Z1, Z2 + Zpara(Z3, Z_ext)) 

    T1 = qubit_radiative_T1(C_q, Zq_env)

    return T1

def qubit_radiative_decay_equivalent_LE_circuit_single_resonator(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50):
    
    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0)

    eff_k = k_readout(C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0, Zline)

    eff_C_ext = get_eff_C_ext_from_k(C1, L1, Zline, eff_k)

    Z1 = Zpara(Zcap(C1, omegas), Zind(L1, omegas)) # 1/(1j*omegas*C1 + 1/(1j*omegas*L1))

    Zg = Zcap(C_g, omegas)

    Z_ext = Zcap(eff_C_ext, omegas) + Zline 

    Zq_env = Zg + Zpara(Z1, Z_ext) 

    T1 = qubit_radiative_T1(C_q, Zq_env)

    return T1

#####################
### get ham terms ###
#####################

def g_coupling(omega_q, C_q, C_g, l_c, l_Gf, l_Gn, phase_vel=3*10**8/2.5, Z0=65):

    ### approx value for qubit readout resonator coupling

    cpw_length = l_c + l_Gf + l_Gn

    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw_length)

    omega_r_val = omega_r(C1, L1)

    val = 0.25 * np.sqrt(omega_q * omega_r_val / (C_q * C1)) * C_g * (omega_q/omega_r_val + omega_r_val/omega_q)

    return val

def J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, receiver_type = 'lambda/4'):

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0, receiver_type = receiver_type)
    J = lumped_model_resonator_coupling(C1, L1, C2, L2, Cg, Lg)

    return J

def k_readout(C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, Zline = 50):

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0)

    J_val = lumped_model_resonator_coupling(C1, L1, C2, L2, Cg, Lg)

    k_ext_val = k_ext_from_LE_model(C2, L2, C_ext, Zline)

    #print('C2 (fF):', C2 * 1e15)
    #print('k_ext_val (MHz):', k_ext_val / (2*np.pi * 1e6))

    omega_1 = omega_r(C1, L1)
    omega_2 = omega_r(C2, L2)

    val = k_readout_from_ham(omega_1, omega_2, k_ext_val, J_val)

    return val

def J_coupling_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, receiver_type = 'lambda/4', simplified = True):
    
    L_G = l_c + l_Gf + l_Gn
    L_P = l_c + l_Rf + l_Rn
    
    omega_r = lambda_quarter_omega(L_G, phase_vel=phase_vel)
    omega_p = lambda_quarter_omega(L_P, phase_vel=phase_vel)

    omega_n = notch_filter_frequency_rule_of_thumb(l_c, l_Gf, l_Rf, phase_vel=phase_vel, Z0=Z0)

    if receiver_type == 'lambda/2':
        omega_n = 2 * omega_n
        omega_p = 2 * omega_p

    C_l = 1/(Z0 * phase_vel)

    if receiver_type == 'lambda/4':

        J_val = J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=phase_vel, Z0=Z0, simplified = simplified)
    
    if receiver_type == 'lambda/2':#
        
        omega_bar = (omega_r + omega_p )/2

        J_val = np.pi **2 / 16 * omega_bar * (omega_r/omega_n - omega_n/omega_r) * (omega_p/omega_n - omega_n/omega_p) * (omega_bar/omega_n - omega_n/omega_bar) * (Cm_per_len /C_l) * np.sin(omega_n * l_c / phase_vel) * 1/(np.cos(omega_n * np.pi / (2*omega_r)) * np.sin(omega_n * np.pi / (omega_p)))

    return J_val

def J_coupling_analytic_by_freqs(omega_r, omega_p, omega_n, l_c, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, simplified = True):

    C_l = 1/(Z0 * phase_vel)

    if simplified:
        ## return the solution ignoring second order terms in the parameter delta_rp/omega_n
        
        omega_bar = (omega_r + omega_p) / 2

        J_val = np.pi **2 / 32 * omega_bar * (omega_bar/omega_n - omega_n/omega_bar)**3 * (Cm_per_len /C_l) * np.sin(omega_n * l_c / phase_vel) * 1/(np.cos(omega_n * np.pi / (2*omega_bar)))**2
    
    else:
        ## return the exact solution to Solgun's simple impedance formula
        J_val = (np.pi * omega_r * omega_p / 8)**2 * (omega_r/omega_n - omega_n/omega_r) * (omega_p/omega_n - omega_n/omega_p) * ((omega_r/omega_n - omega_n/omega_r)/omega_r**3 + (omega_p/omega_n - omega_n/omega_p)/omega_p**3) * (Cm_per_len /C_l) * np.sin(omega_n * l_c / phase_vel) * 1/(np.cos(omega_n * np.pi / (2*omega_r)) * np.cos(omega_n * np.pi / (2*omega_p)))

    return J_val

######################
### User functions ###
######################

def lambda_quarter_omega(cpw_length, phase_vel=3*10**8/2.5):
    # Small utility function to calculate the resonance frequency of an uncoupled lambda/4 resonator

    lambda_line = 4 * cpw_length
    val = 2 * np.pi * phase_vel / lambda_line

    return val

def lumped_elements_J_formula(om1, om2, Z1, Z2, L1, L2):

    return -0.25*np.sqrt(om1*om2/L1/L2)*np.imag(Z1/om1+Z2/om2)

def k_ext_from_LE_model(Cr, Lr, C_ext, Zline = 50):

    k_ext = C_ext**2 * Zline / (Cr**2 * Lr)

    return k_ext

def k_readout_from_ham(omega_1, omega_2, k_ext, J):

    delta = omega_2 - omega_1

    sqrt_term = (k_ext - 2j*delta)**2 - 16*J**2

    val = 0.5 * k_ext * (1 - np.real(np.emath.sqrt(sqrt_term))/k_ext)

    return val

def qubit_radiative_T1(Cq, Zq_env):

    ### using the expression for the radiative decay time found in ref: Houck et al., Phys. Rev. Lett. 101, 080502 (2008)

    Yq_env = 1/Zq_env

    R = 1/ np.abs(np.real(Yq_env))

    T1 = R * Cq

    return T1

def get_eff_Cf_from_J(l_c, l_Gf, l_Gn, l_Rf, l_Rn, J, phase_vel=3*10**8/2.5, Z0=65):

    ### function for getting the value of the coupling capacitance between the resonators, Cf, 
    ### given a particular value of J and assuming that there is no shunt impedance Lf in parallel with Cf

    cpw_length1 = l_c + l_Gf + l_Gn
    cpw_length2 = l_c + l_Rf + l_Rn

    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw_length1)

    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw_length2)

    omega1 = omega_r(C1, L1)
    omega2 = omega_r(C2, L2)

    val = 4 * J * np.sqrt(C1 * C2/(omega1 * omega2)) / (omega1/omega2 + omega2/omega1)

    return val

def get_eff_C_ext_from_k(Cr, Lr, Zline, k_readout):

    val =  np.sqrt(k_readout / (Zline / (Cr**2 * Lr)))
    
    return val

def enhancement_factor_symbolic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, omegas, phase_vel = 3*10**8/2.5):
    ## returns the predicted enhancement to T1 that comes from using the ind-cap coupling design
    ## vs using a direct capacitive coupling design with the same J coupling between the resonators
    ## only valid for omegas near the notch frequency

    omega_r = lambda_quarter_omega(l_c + l_Gf + l_Gn, phase_vel=phase_vel)
    omega_p = lambda_quarter_omega(l_c + l_Rf + l_Rn, phase_vel=phase_vel)
    omega_n = notch_filter_frequency_rule_of_thumb(l_c, l_Gf, l_Rf)

    delta = omegas - omega_n
    omega_bar = (omega_r + omega_p)/ 2

    ## first order
    Z21_ratios =  (2*delta) / (omega_n*((omega_n/omega_bar)**2 -1))

    ## more accurate
    #Z21_ratios =  (omega_n/omegas - omegas/omega_n)/omegas * omega_bar/(omega_n/omega_bar - omega_bar/omega_n)

    T1_enhancement = 1/Z21_ratios**2

    return T1_enhancement

def notch_enhancement_bandwidth_by_omegas(omega_n, omega_r, omega_p, T1_enhancement_fact = 10):

    omega_bar = (omega_r + omega_p)/ 2

    factor = 1/np.sqrt(T1_enhancement_fact)
    
    bandwidth_val = factor * omega_n*(1 - (omega_n/omega_bar)**2)

    return bandwidth_val

def notch_enhancement_bandwidth(l_c, l_Gf, l_Gn, l_Rf, l_Rn, T1_enhancement_fact = 10, phase_vel = 3*10**8/2.5):
    ## returns the predicted bandwidth of the enhancement to T1 that comes from using the ind-cap coupling design
    ## vs using a direct capacitive coupling design with the same J coupling between the resonators

    omega_r = lambda_quarter_omega(l_c + l_Gf + l_Gn, phase_vel=phase_vel)
    omega_p = lambda_quarter_omega(l_c + l_Rf + l_Rn, phase_vel=phase_vel)
    omega_n = notch_filter_frequency_rule_of_thumb(l_c, l_Gf, l_Rf)

    bandwidth_val = notch_enhancement_bandwidth_by_omegas(omega_n, omega_r, omega_p, T1_enhancement_fact = T1_enhancement_fact)

    return bandwidth_val

def find_notch_enhancement_bandwidth(omegas, enhancement_factors, bandwidth_threshold):

    idx = np.argwhere(np.diff(np.sign(enhancement_factors - bandwidth_threshold))).flatten()

    omega_vals = omegas[idx]

    print('omega_vals:', omega_vals)

    if len(omega_vals)<2:
        warnings.warn("Could not find two frequencies that intercept with the bandwidth_threshold so returning a nan")
        enhancement_bandwidth = np.nan
    elif len(omega_vals)>2:
        warnings.warn("Found more than two frequencies that intercept with the bandwidth_threshold so using the first two")
        enhancement_bandwidth = np.abs(omega_vals[1] - omega_vals[0])
    else:
        enhancement_bandwidth = np.abs(omega_vals[1] - omega_vals[0])

    return enhancement_bandwidth