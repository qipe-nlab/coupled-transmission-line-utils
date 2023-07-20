from cmath import phase
import re
import numpy as np
import matplotlib.pyplot as plt
import sys

#######################
### Basic functions ###
#######################

def Zcap(Cval, omega):
	return -1j/(Cval*omega)

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

    fp_exp_end, bp_exp_end = prop_exp_mats(prop_eigs, len)

    mat11 = (Zchar + Z_nb)@T
    mat12 = (Zchar - Z_nb)@T
    mat21 = (Zchar - Z_fb)@T@fp_exp_end
    mat22 = (Zchar + Z_fb)@T@bp_exp_end

    mat = np.block([[mat11, mat12],[mat21, mat22]])

    return mat

def uc_I_sols(defining_matrix, Vs_nb, Vs_fb):

    # solve for the diagonalized currents

    defining_mat_inv = np.linalg.inv(defining_matrix)

    Vs = np.concatenate((Vs_nb, Vs_fb))

    Is = defining_mat_inv@Vs

    Is_split = np.array_split(Is, 2)

    forward_uc_I_sols = Is_split[0]
    backward_uc_I_sols = Is_split[1]

    return forward_uc_I_sols, backward_uc_I_sols 

def V_I_sols(Z_char_mat_sol, T_mat_sol, forward_uc_I_sols, backward_uc_I_sols, prop_eigs, pos):

    # solve for the voltage and current solutions in the original (undiagonalized) frame
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

def voltage_at_source_location_exact(Z0, Zinput_Gn, phase_vel, l_Gn, omega):

    # returns the exact ratio of the input current to the generator line to the voltage at the start of the coupled line section.

    tau_n = l_Gn / phase_vel

    Zinput_n = Zinput(Z0, Zinput_Gn, omega * tau_n)

    I_in = 1

    v_in = Zinput_n * I_in

    v_out, I_out = transmission_line_voltage_current_out(Z0, tau_n, v_in, I_in)

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

    forward_uc_I_sols, backward_uc_I_sols  = uc_I_sols(defining_mat_sol, Vs_nb, Vs_fb)

    V_sols_n, I_sols_n = V_I_sols(Z_char_mat_sol, T_mat_sol, forward_uc_I_sols, backward_uc_I_sols, propagation_eig_sol, 0)

    V_sols_f, I_sols_f = V_I_sols(Z_char_mat_sol, T_mat_sol, forward_uc_I_sols, backward_uc_I_sols, propagation_eig_sol, len)

    return V_sols_n, V_sols_f, I_sols_n, I_sols_f

########################################
### Distributed Element Calculations ###
### 2 coupled transmission lines, WC ###
########################################

### For 2 coupled transmission lines under the assumption of weak coupling
### See: Solution of the Transmission-Line Equations Under the Weak-Coupling Assumption

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

def Z_short_tl(Z0, phase_vel, L, omega):
    tau = L/phase_vel
    return 1j*Z0*np.tan(tau * omega)

def Z_open_tl(Z0, phase_vel, L, omega):
    tau = L/phase_vel
    return -1j*Z0/np.tan(tau * omega)

def voltage_transmission_coupled_lines(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega):

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

    Z_Gn = 0 # 1j * Z0 * (np.tan(tau_G*omega)*np.cos(tau_Gn*omega) - np.sin(tau_Gn*omega))  / (np.tan(tau_G*omega)*np.sin(tau_Gn*omega) + np.cos(tau_Gn*omega)) 

    gamma_Rn = reflection_ceofficient(Z_Rn, Z0)

    gamma_Rf = reflection_ceofficient(Z_Rf, Z0)

    gamma_Gn = reflection_ceofficient(Z_Gn, Z0)

    gamma_Gf = reflection_ceofficient(Z_Gf, Z0)

    F_G = F_val(gamma_Gn, gamma_Gf, tau_c, omega)

    F_R = F_val(gamma_Rn, gamma_Rf, tau_c, omega)

    A_value = A_val(K_plus, K_minus, tau_c, gamma_Gf, omega)

    B_value = B_val(K_plus, K_minus, tau_c, gamma_Gf, omega)

    val = (Z_Rn/(Z_Rn + Z0_c)) * (A_value + B_value * gamma_Rf * np.exp(-1j*omega*tau_c)) * F_G * F_R * (Z0_c/(Z0_c + Z_Gn))
    
    #val = (A_value + B_value * gamma_Rf * np.exp(1j*omega*tau_c))

    #val = (A_value + B_value * gamma_Rf * np.exp(1j*omega*tau_c)) *  F_G * F_R

    return val

def Z_trans_along_shorted_tl(Z0, phase_vel, L, z, omega):

    tau = L/phase_vel

    val = 1j * Z0 * (np.tan(tau*omega)*np.cos(tau*z/L * omega) - np.sin(tau*z/L * omega)) 

    return val

def voltage_at_source_location(Z0, phase_vel, Cm_per_len, l_c, l_Gf, l_Gn, omega):

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

    v_out, I_out = transmission_line_voltage_current_out(Z0, tau_n, v_in, I_in)

    val = v_out

    return val

def find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=3*10**8/2.5, Z0=65, min_search=3*2*np.pi*10**9, max_search=12*2*np.pi*10**9, search_spacing=(0.01*2*np.pi*10**9)):
        
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

def lumped_model_C_and_L(phase_vel, Z0, cpw__length):

    tl_C_val, tl_L_val = transmission_line_C_and_L(phase_vel, Z0)

    C_val = tl_C_val * cpw__length / 2
    L_val = 8 * tl_L_val * cpw__length / np.pi**2 

    return C_val, L_val

def lumped_model_Cg_and_Lg(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len):

    omega_f = find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0, 1e9*2*np.pi, 10e9*2*np.pi, 5e4*2*np.pi)
    Z0_f = find_notch_filter_char_impedance(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega_f, phase_vel=phase_vel, Z0=Z0)

    Cg = 1/(omega_f*Z0_f)
    Lg =  Z0_f / omega_f

    return Cg, Lg

def get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65):
        
    cpw__length1 = l_c + l_Gf + l_Gn
    cpw__length2 = l_c + l_Rf + l_Rn
    
    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1)
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2)
    Cg, Lg = lumped_model_Cg_and_Lg(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)
    
    Cg_basic = Cm_per_len *l_c
    Lg_basic = L1*L2/ (Lm_per_len * l_c)

    return C1, L1, C2, L2, Cg, Lg

def lumped_model_Z_transmission(omega, C1, L1, C2, L2, Cg, Lg):

    Z_res1 = Zres(C1, L1, omega)
    Z_res_coupler = Zres(Cg, Lg, omega)
    Z_res2 = Zres(C2, L2, omega)

    val = Z_res1 * Z_res2 / Z_res_coupler * 1/(1 + (Z_res1 + Z_res2)/Z_res_coupler)

    return val 

def lumped_model_resonator_coupling(C1, L1, C2, L2, Cg, Lg):

    omega_1 = omega_r(C1, L1)
    omega_2 = omega_r(C2, L2)

    Z21_omega_1 = lumped_model_Z_transmission(omega_1, C1, 1, C2, 1, Cg, Lg)
    Z21_omega_1 = lumped_model_Z_transmission(omega_2, C1, 1, C2, 1, Cg, Lg)

    val = -0.25*(omega_1**3*omega_2**3*C1*C2)**0.5*(np.imag(Z21_omega_1)/omega_1 + np.imag(Z21_omega_1)/omega_2)

    return val

def find_notch_filter_char_impedance(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f, phase_vel=3*10**8/2.5, Z0=65):

    cpw__length1 = l_c + l_Gf + l_Gn
    cpw__length2 = l_c + l_Rf + l_Rn

    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1)
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2)

    Zres_1_omega_f = Zres(C1, L1, omega_f)
    Zres_2_omega_f = Zres(C2, L2, omega_f)

    delta_omega = 10 * 2*np.pi # 10 Hz

    Z_transfer_plus = Z_transfer_weak_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f + delta_omega/2, phase_vel=phase_vel, Z0=Z0)
    Z_transfer_minus = Z_transfer_weak_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f - delta_omega/2, phase_vel=phase_vel, Z0=Z0)

    grad_Ztransfer = (Z_transfer_plus - Z_transfer_minus)/(delta_omega)

    val = np.real(1j*2*Zres_1_omega_f * Zres_2_omega_f / (omega_f * grad_Ztransfer))

    return val

###################################
### Z matrix transfer functions ###
###################################

def Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=3*10**8/2.5, Z0=65):

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
    Rn_tl_voltage_scale_factor_arr = []
    Z_input_exact_arr = []

    if len(omegas) == 0:
        omegas = [omegas]

    for omega in omegas:

        Z_nb = np.diag(np.array([0, Z_open_tl(Z0, phase_vel, l_Rn, omega)]))
        Z_fb = np.diag(np.array([Z_short_tl(Z0, phase_vel, l_Gf, omega), Z_short_tl(Z0, phase_vel, l_Rf, omega)]))

        Vs_nb = np.array([1,0])
        Vs_fb = np.array([0,0])

        V_sols_n, V_sols_f, I_sols_n, I_sols_f = V_I_solutions_sym_three(C_Maxwell_test, L_test, Z_nb, Z_fb, Vs_nb, Vs_fb, l_c, omega)

        V_sols_n_arr.append(V_sols_n)
        V_sols_f_arr.append(V_sols_f)

        I_sols_n_arr.append(I_sols_n)
        I_sols_f_arr.append(I_sols_f)

        Z_input_exact = V_sols_n[0] / I_sols_n[0]

        Z_input_exact_arr.append(Z_input_exact)

        #print('Z_input_exact:', Z_input_exact)

        V_from_input_I = voltage_at_source_location_exact(Z0, Z_input_exact, phase_vel, l_Gn, omega)

        #print('V_from_input_I:', V_from_input_I)

        V_from_input_I_arr.append(V_from_input_I)

        tau_Rn = l_Rn/phase_vel

        Rn_tl_voltage_scale_factor = 1/np.cos(tau_Rn*omega)

        Rn_tl_voltage_scale_factor_arr.append(Rn_tl_voltage_scale_factor)

    # print('V_sols_n_arr:', V_sols_n_arr)
    # print('V_sols_f_arr:', V_sols_f_arr)

    V_sols_n_arr = np.array(V_sols_n_arr)
    V_sols_f_arr = np.array(V_sols_f_arr)
    V_from_input_I_arr = np.array(V_from_input_I_arr)
    Rn_tl_voltage_scale_factor_arr = np.array(Rn_tl_voltage_scale_factor_arr)

    Z_input_exact_arr = np.array(Z_input_exact_arr)

    ## take the imag part to remove tiny numerical errors

    Z_transfer_exact = np.imag(Rn_tl_voltage_scale_factor_arr * V_from_input_I_arr * V_sols_n_arr[:, -1])

    return Z_transfer_exact

def Z_transfer_weak_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega, phase_vel=3*10**8/2.5, Z0=65):

    L_G = l_Gn + l_c + l_Gf

    #L_G = l_c + l_Gf

    #Z_trans_along_shorted_tr_val = Z_trans_along_shorted_tl(Z0, phase_vel, L_G, l_Gn, omega)

    Z_trans_along_shorted_tr_val = voltage_at_source_location(Z0, phase_vel, Cm_per_len, l_c, l_Gf, l_Gn, omega)

    voltage_transmission_coupled_lines_val = voltage_transmission_coupled_lines(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega)

    tau_Rn = l_Rn/phase_vel

    Rn_tl_voltage_scale_factor = 1/np.cos(tau_Rn*omega)

    val = Z_trans_along_shorted_tr_val * voltage_transmission_coupled_lines_val * Rn_tl_voltage_scale_factor

    return val

def Z_transfer_equivalent_LE_circuit(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega, phase_vel=3*10**8/2.5, Z0=65):

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0)

    val = lumped_model_Z_transmission(omega, C1, L1, C2, L2, Cg, Lg)

    return val

#######################################
### qubit radiative decay functions ###
#######################################

def qubit_radiative_decay_sym_3_lines_exact(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50):

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
    Z_input_exact_arr = []

    if len(omegas) == 0:
        omegas = [omegas]

    for omega in omegas:

        Z_env = Zline + Zcap(C_ext, omega)

        Z_nb = np.diag(np.array([0, Z_input_tl(Z0, Z_env, phase_vel, l_Rn, omega)]))
        Z_fb = np.diag(np.array([Z_short_tl(Z0, phase_vel, l_Gf, omega), Z_short_tl(Z0, phase_vel, l_Rf, omega)]))

        Vs_nb = np.array([1,0])
        Vs_fb = np.array([0,0])

        V_sols_n, V_sols_f, I_sols_n, I_sols_f = V_I_solutions_sym_three(C_Maxwell_test, L_test, Z_nb, Z_fb, Vs_nb, Vs_fb, l_c, omega)

        V_sols_n_arr.append(V_sols_n)
        V_sols_f_arr.append(V_sols_f)

        I_sols_n_arr.append(I_sols_n)
        I_sols_f_arr.append(I_sols_f)

        Z_input_exact = V_sols_n[0] / I_sols_n[0]

        Z_input_exact_arr.append(Z_input_exact)

    Z_input_total_exact = np.array([Z_input_tl(Z0, Z_input_exact_arr[i], phase_vel, l_Gn, omega) for i, omega in enumerate(omegas)])

    Zq_env = Zcap(C_g, omegas) + Z_input_total_exact

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

def J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65):

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0)
    J = lumped_model_resonator_coupling(C1, L1, C2, L2, Cg, Lg)

    return J

def k_readout(C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel=3*10**8/2.5, Z0=65, Zline = 50):

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0)

    J_val = lumped_model_resonator_coupling(C1, L1, C2, L2, Cg, Lg)

    k_ext_val = k_ext_from_LE_model(C2, L2, C_ext, Zline)

    print('C2 (fF):', C2 * 1e15)
    print('k_ext_val (MHz):', k_ext_val / (2*np.pi * 1e6))

    omega_1 = omega_r(C1, L1)
    omega_2 = omega_r(C2, L2)

    val = k_readout_from_ham(omega_1, omega_2, k_ext_val, J_val)

    return val

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
