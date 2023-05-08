import numpy as np
import matplotlib.pyplot as plt
import sys

import common_formulas as cf



########################################
### Distributed Element Calculations ###
###     General Numerical Solver     ###
########################################

# General exact solver for n coupled lines
# Only assumption is that lines are lossless
# See: Analysis of Multiconductor Transmission Lines, 2E, Chapter 7

# NOTE: This approach is not used for the dimensioning algorithm
#       See r_r_formulas.py instead

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

    sqrt_U_eig = np.sqrt(U_eig)

    sqrt_U_eig_mat = np.diag(sqrt_U_eig)

    M_mat = sqrt_U_eig_mat @ (U_mat.T) @ L_mat @ U_mat @ sqrt_U_eig_mat

    eig, mat = np.linalg.eigh(M_mat)

    return mat, eig

def T_mat(U_mat, S_mat, U_eig):

    # for diagonalizing CL matrix

    sqrt_U_eig = np.sqrt(U_eig)

    mat = U_mat@np.diag(sqrt_U_eig)@S_mat

    return mat

def T_v_mat(U_mat, S_mat, U_eig):

    # for diagonalizing CL matrix

    sqrt_U_eig = np.sqrt(U_eig)

    mat = U_mat@np.diag(1/sqrt_U_eig)@S_mat

    return mat

def T_inv_mat(U_mat, S_mat, U_eig):

    mat = (T_v_mat(U_mat, S_mat, U_eig)).T

    return mat 

def propagation_eig(S_eig, omega):

    # diagonalized YZ matrix elements

    eig = 1j*omega*np.sqrt(S_eig)

    return eig

def Z_char_mat(U_mat, S_mat, U_eig, S_eig):

    sqrt_U_eig = np.sqrt(U_eig)

    sqrt_S_eig = np.sqrt(S_eig)

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

    Z_char_plus = np.sqrt((Z_mat[0,0] + Z_mat[0,1])/(Y_mat[0,0] + Y_mat[0,1]))
    Z_char_minus =  np.sqrt((Z_mat[0,0] - Z_mat[0,1])/(Y_mat[0,0] - Y_mat[0,1]))

    mat = 0.5 * np.array([[Z_char_plus + Z_char_minus, Z_char_plus - Z_char_minus],[Z_char_plus - Z_char_minus, Z_char_plus + Z_char_minus]])

    return mat

def propagation_eig_sym_three(Z_mat, Y_mat):

    # diagonalized YZ matrix elements

    eig1 = np.sqrt((Z_mat[0,0] + Z_mat[0,1])*(Y_mat[0,0] + Y_mat[0,1]))
    
    eig2 = np.sqrt((Z_mat[0,0] - Z_mat[0,1])*(Y_mat[0,0] - Y_mat[0,1]))

    eig = np.array([eig1, eig2])

    return eig

# basic util functions

def Z_short_tl(Z0, phase_vel, L, omega):

    # impedance of a shorted TL section 

    tau = L/phase_vel

    val = 1j*Z0*np.tan(tau * omega)

9    return val 

def Z_open_tl(Z0, phase_vel, L, omega):

    # impedance of an open TL section 

    tau = L/phase_vel

    val = -1j*Z0/np.tan(tau * omega)

    return val

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

    Z_char_mat_sol = Z_char_mat(U_mat_sol, S_mat_sol, U_eig_sol, S_eig_sol, L_mat)

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
