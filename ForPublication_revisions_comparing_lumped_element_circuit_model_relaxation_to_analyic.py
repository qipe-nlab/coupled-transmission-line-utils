import numpy as np
import matplotlib.pyplot as plt
import json
import sys
from scipy.optimize import minimize

#### functions

### effective kappa_ro

def hybridized_eigenvalues(freq_r, freq_p, kappa_p , J, chi, state = 'g'):

    filter_external_decay_rate = kappa_p
    J_coupling = J
    chi_r = chi
    dressed_cavity_frequency = freq_r
    dressed_filter_frequency = freq_p

    if state == 'g':
        dressed_cavity_freq = dressed_cavity_frequency

    elif state == 'e':
        dressed_cavity_freq = dressed_cavity_frequency + 2*chi_r
    
    delta_rp = dressed_cavity_freq - dressed_filter_frequency

    eigenmode_common = (dressed_cavity_freq + dressed_filter_frequency - 1j*filter_external_decay_rate/2)/2
    eigenmode_diff = 0.5*np.emath.sqrt((delta_rp + 1j*filter_external_decay_rate/2)**2 + 4*J_coupling**2)

    low_mode = eigenmode_common - eigenmode_diff
    high_mode = eigenmode_common + eigenmode_diff

    return low_mode, high_mode

def hybridized_modes_kappas(freq_r, freq_p, kappa_p, J, chi, state = 'g'):

    low_mode, high_mode = hybridized_eigenvalues(freq_r, freq_p, kappa_p, J, chi, state = state)

    kappa_low = -2*np.imag(low_mode)
    kappa_high = -2*np.imag(high_mode)

    return kappa_low, kappa_high

### 

def Zpara(Z1, Z2):
    return 1/(1/Z1+1/Z2)

def Zcap(Cval, omega):

    if Cval == 0:
        val = 1j*1e30
    else:
        val = -1j/(Cval*omega)

    return val

def Zind(Lval, omega):
    return 1j*(Lval*omega)

def transmission_line_C_and_L(phase_vel, Z0):

    C_val = 1/(phase_vel*Z0)
    L_val = Z0 / phase_vel

    return C_val, L_val

def lumped_model_C_and_L(phase_vel, Z0, cpw_length, res_type = 'lambda/4'):

    tl_C_val, tl_L_val = transmission_line_C_and_L(phase_vel, Z0)

    if res_type == 'lambda/4':
        C_val = tl_C_val * cpw_length / 2
        L_val = 8 * tl_L_val * cpw_length / np.pi**2 

    elif res_type == 'lambda/2':
        L_val = tl_L_val * cpw_length / 2
        C_val = 2*tl_C_val * cpw_length / np.pi**2 

    return C_val, L_val

### initial cicrcuit param functions

def lumped_model_C_and_L_from_freq(phase_vel, Z0, omega_res, res_type = 'lambda/4'):

    length = np.pi*phase_vel/(2*omega_res)
    
    #print('length:', length)

    C_val, L_val = lumped_model_C_and_L(phase_vel, Z0, length, res_type = res_type)

    return C_val, L_val

def get_eff_Cg_from_g(omega_q, omega_r, C_q, Cr, g):

    ### approx value for qubit readout resonator coupling capacitance
    
    val = 4 * g/((omega_q/omega_r + omega_r/omega_q) * np.sqrt(omega_q * omega_r / (C_q * Cr)))
    
    return val

def get_eff_Cf_from_J(C1, L1, C2, L2, omega1, omega2, J):

    ### function for getting the value of the coupling capacitance between the resonators, Cf, 
    ### given a particular value of J and assuming that there is no shunt impedance Lf in parallel with Cf

    val = 4 * J * np.sqrt(C1 * C2/(omega1 * omega2)) / (omega1/omega2 + omega2/omega1)

    return val

def get_eff_C_ext_from_k(Cr, Lr, Zline, k_readout):

    val =  np.sqrt(k_readout / (Zline / (Cr**2 * Lr)))
    
    return val

### exact cicrcuit param functions

def generate_bounds(args):
    
    bounds = [(arg/2, arg*2) for arg in args]
    
    return bounds

def omega_r_func(C_r, C_g, L_r):
    
    val = 1/np.sqrt((C_r + C_g)*L_r)
    
    return val

def omega_p_func(C_p, C_k, L_p):
    
    val = 1/np.sqrt((C_p + C_k)*L_p)
    
    return val

def g_func(C_q, C_r, C_g, omega_q, omega_r):

    val = 0.5 * C_g/np.sqrt((C_r + C_g)*(C_q + C_g)) * np.sqrt(omega_q*omega_r)
    
    return val
    
def J_cap(C_r, C_p, C_g, C_J, C_k, omega_r, omega_p):
    
    val = 0.5 * C_J/np.sqrt((C_r + C_g)*(C_p + C_k)) * np.sqrt(omega_r*omega_p)
    
    return val

def J_MTL(C_r, C_p, C_g, L_J, C_J, C_k, omega_r, omega_p, Z0):
    
    omega_n = 1/np.sqrt(L_J*C_J)
    
    Z0n = np.sqrt(L_J/C_J)
    
    val = 2 * Z0 / (Z0n * np.pi) * (omega_r*omega_p)**0.5 * ((omega_r*omega_p)**0.5/omega_n - omega_n/(omega_r*omega_p)**0.5)
    
    return val

def kappa_p_func(C_p, L_p, C_k, Z_env):
    
    val = Z_env/L_p * (C_k/(C_k + C_p))**2
    
    return val

def mini_func_with_notch(x, omega_q, omega_r, omega_p, omega_n, g, J, k, C_q, Z0, Z_env):
    
    L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k = tuple(x)
    
    ### rescale variables to aid fitter
    L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k, C_q = tuple(np.array([L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k, C_q])*1e10)
    omega_q, omega_r, omega_p, omega_n, g, J, k = tuple(np.array([omega_q, omega_r, omega_p, omega_n, g, J, k])*1e-10)
    
    omega_r_resi = np.abs(omega_r_func(C_r, C_g, L_r) - omega_r)**2
    omega_p_resi = np.abs(omega_p_func(C_p, C_k, L_p) - omega_p)**2
    g_resi = np.abs(g_func(C_q, C_r, C_g, omega_q, omega_r) - g)**2
    J_MTL_resi =  np.abs(J_MTL(C_r, C_p, C_g, L_J, C_J, C_k, omega_r, omega_p, Z0) - J)**2
    kappa_p_resi =  np.abs(kappa_p_func(C_p, L_p, C_k, Z_env) - k)**2
    
    resi = omega_r_resi + omega_p_resi + J_MTL_resi + kappa_p_resi
    
    resi = resi*1e6
        
    return resi

def mini_func_without_notch(x, omega_q, omega_r, omega_p, g, J, k, C_q, Z0, Z_env):
    
    L_r, L_p, C_r, C_p, C_g, C_J, C_k = tuple(x)
    
    ### rescale variables to aid fitter
    L_r, L_p, C_r, C_p, C_g, C_J, C_k, C_q = tuple(np.array([L_r, L_p, C_r, C_p, C_g, C_J, C_k, C_q])*1e10)
    omega_q, omega_r, omega_p, g, J, k = tuple(np.array([omega_q, omega_r, omega_p, g, J, k])*1e-10)

    omega_r_resi = np.abs(omega_r_func(C_r, C_g, L_r) - omega_r)**2
    omega_p_resi = np.abs(omega_p_func(C_p, C_k, L_p) - omega_p)**2
    g_resi = np.abs(g_func(C_q, C_r, C_g, omega_q, omega_r) - g)**2
    J_cap_resi =  np.abs(J_cap(C_r, C_p, C_g, C_J, C_k, omega_r, omega_p) - J)**2
    kappa_p_resi =  np.abs(kappa_p_func(C_p, L_p, C_k, Z_env) - k)**2
    
    resi = omega_r_resi + omega_p_resi + J_cap_resi + kappa_p_resi
    
    resi = resi*1e6
    
    return resi

def mini_func_without_filter(x, omega_q, omega_r, g, k, C_q, Z0, Z_env):
    
    L_r, C_r, C_g, C_k = tuple(x)

    ### rescale variables to aid fitter
    L_r, C_r, C_g, C_k, C_q = tuple(np.array([L_r, C_r, C_g, C_k, C_q])*1e10)
    omega_q, omega_r, g, k = tuple(np.array([omega_q, omega_r, g, k])*1e-10)
    
    omega_r_resi = np.abs(omega_r_func(C_r, C_g, L_r) - omega_r)**2
    g_resi = np.abs(g_func(C_q, C_r, C_g, omega_q, omega_r) - g)**2
    kappa_ro_resi =  np.abs(kappa_p_func(C_r, L_r, C_k, Z_env) - k)**2

    resi = omega_r_resi + g_resi + kappa_ro_resi

    resi = resi*1e6
    
    return resi

def get_initial_circuit_params_with_notch(omega_q, omega_r, omega_p, omega_n, g, J, k, C_q, phase_vel, Z0, Z_env):
    
    C_r, L_r = lumped_model_C_and_L_from_freq(phase_vel, Z0, omega_r)
    C_p, L_p = lumped_model_C_and_L_from_freq(phase_vel, Z0, omega_p)
    C_g = get_eff_Cg_from_g(omega_q, omega_r, C_q, C_r, g)
    C_k = get_eff_C_ext_from_k(C_r, L_r, Z_env, k)

    Z0n = 2 * Z0 / (J * np.pi) * (omega_r*omega_p)**0.5 * ((omega_r*omega_p)**0.5/omega_n - omega_n/(omega_r*omega_p)**0.5)

    C_J  = 1/(Z0n * omega_n)
    L_J = Z0n / omega_n

    return L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k

def get_initial_circuit_params_without_notch(omega_q, omega_r, omega_p, g, J, k, C_q, phase_vel, Z0, Z_env):
    
    C_r, L_r = lumped_model_C_and_L_from_freq(phase_vel, Z0, omega_r)
    C_p, L_p = lumped_model_C_and_L_from_freq(phase_vel, Z0, omega_p)  
    C_g = get_eff_Cg_from_g(omega_q, omega_r, C_q, C_r, g)
    C_k = get_eff_C_ext_from_k(C_r, L_r, Z_env, k)
    C_J = get_eff_Cf_from_J(C_r, L_r, C_p, L_p, omega_r, omega_p, J)
    
    return L_r, L_p, C_r, C_p, C_g, C_J, C_k

def get_initial_circuit_params_without_filter(omega_q, omega_r, g, k, C_q, phase_vel, Z0, Z_env):
    
    C_r, L_r = lumped_model_C_and_L_from_freq(phase_vel, Z0, omega_r)
    C_g = get_eff_Cg_from_g(omega_q, omega_r, C_q, C_r, g)
    C_k = get_eff_C_ext_from_k(C_r, L_r, Z_env, k)
    
    return L_r, C_r, C_g, C_k

def find_circuit_params_with_notch(omega_q, omega_r, omega_p, omega_n, g, J, k, C_q, phase_vel, Z0, Z_env):
    
    L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k = get_initial_circuit_params_with_notch(omega_q, omega_r, omega_p, omega_n, g, J, k, C_q, phase_vel, Z0, Z_env)

    x0 = np.array([L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k])
    
    bound_vals = generate_bounds(x0)
    
    #print('bound_vals:', bound_vals)
    
    res = minimize(mini_func_with_notch, x0, args = (omega_q, omega_r, omega_p, omega_n, g, J, k, C_q, Z0, Z_env), bounds = bound_vals, method = 'Nelder-Mead', tol = 1e9)

    #print('res:', res)
    
    params = res.x
    
    #print('params:', params)

    L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k = tuple(params)
        
    return L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k

def find_circuit_params_without_notch(omega_q, omega_r, omega_p, g, J, k, C_q, phase_vel, Z0, Z_env):
    
    ## initial circuit params
    L_r, L_p, C_r, C_p, C_g, C_J, C_k = get_initial_circuit_params_without_notch(omega_q, omega_r, omega_p, g, J, k, C_q, phase_vel, Z0, Z_env)

    x0 = np.array([L_r, L_p, C_r, C_p, C_g, C_J, C_k])
    
    bound_vals = generate_bounds(x0)
    
    res = minimize(mini_func_without_notch, x0, args = (omega_q, omega_r, omega_p, g, J, k, C_q, Z0, Z_env), bounds = bound_vals, method = 'Nelder-Mead', tol = 1e9)
    
    #print('res:', res)
    
    params = res.x
    
    #print('params:', params)

    L_r, L_p, C_r, C_p, C_g, C_J, C_k = tuple(params)
    
    return L_r, L_p, C_r, C_p, C_g, C_J, C_k

def find_circuit_params_without_filter(omega_q, omega_r, g, k, C_q, phase_vel, Z0, Z_env):
    
    ## initial circuit params
    L_r, C_r, C_g, C_k = get_initial_circuit_params_without_filter(omega_q, omega_r, g, k, C_q, phase_vel, Z0, Z_env)

    x0 = np.array([L_r, C_r, C_g, C_k])
    
    bound_vals = generate_bounds(x0)
    
    res = minimize(mini_func_without_filter, x0, args = (omega_q, omega_r, g, k, C_q, Z0, Z_env), bounds = bound_vals, method = 'Nelder-Mead', tol = 1e9)
        
    params = res.x
    
    #print('params:', params)

    L_r, C_r, C_g, C_k = tuple(params)
    
    return L_r, C_r, C_g, C_k

### T1 limit

def qubit_radiative_T1(Cq, Zq_env):

    ### using the expression for the radiative decay time found in ref: Houck et al., Phys. Rev. Lett. 101, 080502 (2008)

    Yq_env = 1/Zq_env

    R = 1/ np.abs(np.real(Yq_env))

    T1 = R * Cq

    return T1

def T1_limit_without_notch(omegas, L_r, L_p, C_r, C_p, C_g, C_J, C_k, with_shunt = True, paper_form = True):
    
    #omegas = np.linspace(freqs_sim[0], freqs_sim[-1], 250)* 2 * np.pi*1e9

    Z1 = 1/(1j*omegas*C_r + 1/(1j*omegas*L_r))
    Z2 = 1/(1j*omegas*C_J)
    Z3 = 1/(1j*omegas*C_p + 1/(1j*omegas*L_p))

    Zg = Zcap(C_g, omegas)

    C_shunt = 230e-15 #+ 20*3e-15 # testing adding the effective shunt from other resonators
    L_shunt = 1.01e-9

    Zline = 50
    
    Zshunt = Zpara(Zcap(C_shunt, omegas), Zind(L_shunt, omegas))

    if paper_form:
        
        Z11 = Zpara(Z1, Z2 + Z3)
        Z22 = Zpara(Z3, Z2 + Z1)
        Z21 = Z11/(Z2 + Z3) * Z3
        
        if with_shunt:
            Z_ext = Zcap(C_k, omegas) + Zshunt/(1+np.abs(Zshunt/Zline)**2)
            Zline = Zline/(1+np.abs(Zline/Zshunt)**2)    
        else:
            Z_ext = Zcap(C_k, omegas)
        
        Re_Yin = Zline * np.abs(Z21)**2 / np.abs((Z11 + Zg)*(Z22 + Z_ext + Zline))**2
        T1_without_notch = Cq / Re_Yin
    
    else:
        if with_shunt:
            Z_ext = Zcap(C_k, omegas) + Zpara(Zshunt, Zline)
        else:
            Z_ext = Zcap(C_k, omegas) + Zline
    
        Zq_env = Zg + Zpara(Z1, Z2 + Zpara(Z3, Z_ext)) 

        T1_without_notch = qubit_radiative_T1(Cq, Zq_env)

    return T1_without_notch

def T1_limit_with_notch(omegas, L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k, with_shunt = True, paper_form = True):
    
    #omegas = np.linspace(freqs_sim[0], freqs_sim[-1], 250)* 2 * np.pi*1e9

    #print('notch_freq:', 1/np.sqrt(L_J*C_J) / (2*np.pi*1e6))
    
    Z1 = 1/(1j*omegas*C_r + 1/(1j*omegas*L_r))
    Z2_with_notch = Zpara(Zcap(C_J, omegas), Zind(L_J, omegas))
    Z3 = 1/(1j*omegas*C_p + 1/(1j*omegas*L_p))
    Zg = Zcap(C_g, omegas)

    Zline = 50
    
    #Z_ext = Zcap(C_k, omegas)  + Zline
    
    C_shunt = 230e-15 #+ 20*3e-15 # testing adding the effective shunt from other resonators
    L_shunt = 1.01e-9
    Zshunt = Zpara(Zcap(C_shunt, omegas), Zind(L_shunt, omegas))        

    if paper_form:
        
        Z11 = Zpara(Z1, Z2_with_notch + Z3)
        Z22 = Zpara(Z3, Z2_with_notch + Z1)
        Z21 = Z11/(Z2_with_notch + Z3) * Z3
        
        if with_shunt:
            Z_ext = Zcap(C_k, omegas) + Zshunt/(1+np.abs(Zshunt/Zline)**2)
            Zline = Zline/(1+np.abs(Zline/Zshunt)**2)    
        else:
            Z_ext = Zcap(C_k, omegas)
        
        Re_Yin = Zline * np.abs(Z21)**2 / np.abs((Z11 + Zg)*(Z22 + Z_ext + Zline))**2
        T1_with_notch = Cq / Re_Yin
        
    else:

        if with_shunt:
            Z_ext = Zcap(C_k, omegas) + Zpara(Zshunt, Zline)
        else:
            Z_ext = Zcap(C_k, omegas) + Zline

        Zq_env_with_notch = Zg + Zpara(Z1, Z2_with_notch + Zpara(Z3, Z_ext))

        T1_with_notch = qubit_radiative_T1(Cq, Zq_env_with_notch)

    return T1_with_notch

def T1_limit_without_filter(omegas, L_r, C_r, C_g, C_k, with_shunt = True, paper_form = True):
    
    Zres = 1/(1j*omegas*C_r + 1/(1j*omegas*L_r))
    Zg = Zcap(C_g, omegas)

    C_shunt = 230e-15 #+ 20*3e-15 # testing adding the effective shunt from other resonators
    L_shunt = 1.01e-9

    Zline = 50
    
    Zshunt = Zpara(Zcap(C_shunt, omegas), Zind(L_shunt, omegas))

    if paper_form:
        
        if with_shunt:
            Z_ext = Zcap(C_k, omegas) + Zshunt/(1+np.abs(Zshunt/Zline)**2)
            Zline = Zline/(1+np.abs(Zline/Zshunt)**2) 
        else:
            Z_ext = Zcap(C_k, omegas)
        
        Re_Yin = np.real(1/(Zg + Zpara(Zres, Z_ext + Zline)))
        T1_without_filter = Cq / Re_Yin
    
    else:
        raise ValueError('paper_form = False not implemented yet')
        # if with_shunt:
        #     Z_ext = Zcap(C_k, omegas) + Zpara(Zshunt, Zline)
        # else:
        #     Z_ext = Zcap(C_k, omegas) + Zline
    
        # Zq_env = Zg + Z1

        # T1_without_filter = qubit_radiative_T1(Cq, Zq_env)

    return T1_without_filter

def T1_enhancement_from_filter_factor(omega_q, omega_r, omega_p,  J, kappa_p, kappa_ro):
    
    omega_rp = (omega_r + omega_p)/2
        
    val = 0.25 * (omega_rp/J)**2 * ((omega_rp/omega_q)**2 - 1)**2 * (kappa_ro/kappa_p)
    
    return val

def T1_enhancement_factor(omega_q, omega_r, omega_p, omega_n, rounded = False):

    ## to MHz
    omega_q = omega_q/(2*np.pi*1e6)
    
    omega_r, omega_p, omega_n = tuple(np.array([omega_r, omega_p, omega_n ])/(2*np.pi*1e6))
    
    omega_rp = (omega_r + omega_p)/2
        
    delta_qn = omega_q - omega_n
    
    if rounded:
        # round to nearest 10 MHz
        delta_qn = np.round(delta_qn,-1)
    
    val = 1/4 * omega_q**2/delta_qn**2 * (1 - (omega_n/omega_rp)**2)**2 
    
    return val

def T1_total_enhancement_factor(omega_q, omega_r, omega_p, omega_n, J, kappa_p, kappa_ro):
    
    T1_enhancement_val = T1_enhancement_factor(omega_q, omega_r, omega_p, omega_n)
    
    T1_enhancement_from_filter_val = T1_enhancement_from_filter_factor(omega_q, omega_r, omega_p, J, kappa_p, kappa_ro)
    
    val = T1_enhancement_val * T1_enhancement_from_filter_val
    
    return val

qubit = 'Q6'

if qubit == 'Q4':
    Cq = 63.8e-15
    freq_q = 8189.48 * 2 * np.pi*1e6 
    freq_r = 10386 * 1e6 * 2*np.pi
    freq_p = 10407 * 1e6 * 2*np.pi
    freq_n = 8.01 * 2 * np.pi * 1e9
    kappa_p = 81.4 * 1e6 * 2*np.pi # 95.42478306 * 1e6 * 2*np.pi #
    kappa_ro = 34 * 2*np.pi*1e6
    J = 39.1 * 1e6 * 2*np.pi
    chi = -9.9 * 1e6 * 2*np.pi
    g = 423 * 2 * np.pi*1e6 

if qubit == 'Q5':
    Cq = 50.6e-15
    freq_q = 8979.8 * 2 * np.pi * 1e6
    freq_r = 10665 * 1e6 * 2*np.pi
    freq_p = 10704 * 1e6 * 2*np.pi
    freq_n = 9.32 * 2 * np.pi * 1e9
    kappa_p = 93.5 * 1e6 * 2*np.pi   # 81.43577923 * 1e6 * 2*np.pi #
    kappa_ro = 19 * 2*np.pi*1e6
    J = 26.2 * 1e6 * 2*np.pi
    chi = -8.3  * 1e6 * 2*np.pi
    g = 275 * 2 * np.pi*1e6
    
if qubit == 'Q6':
    Cq = 51.9e-15
    freq_q = 9045.61 * 2 * np.pi*1e6 
    freq_r = 10540 * 1e6 * 2*np.pi
    freq_p = 10566 * 1e6 * 2*np.pi
    freq_n = 8.775 * 2 * np.pi * 1e9 # 8.91 * 2 * np.pi * 1e9 #
    kappa_p = 66.7  * 1e6 * 2*np.pi # 81.34889282 * 1e6 * 2*np.pi #
    kappa_ro = 24 * 2*np.pi*1e6
    J = 30.9 * 1e6 * 2*np.pi
    chi = -10.5  * 1e6 * 2*np.pi    
    g = 280 * 2 * np.pi*1e6
    
if qubit == 'Q7':
    Cq = 65.1e-15
    freq_q = 8032.04 * 2 * np.pi*1e6 
    freq_r = 10250 * 1e6 * 2*np.pi
    freq_p = 10232 * 1e6 * 2*np.pi
    freq_n = 8.006 * 2 * np.pi * 1e9 #8.005 * 2 * np.pi * 1e9
    kappa_p = 97.6 * 1e6 * 2*np.pi # 85.59018588 * 1e6 * 2*np.pi #
    kappa_ro = 42 * 2*np.pi*1e6
    J = 36.1 * 1e6 * 2*np.pi
    chi = -9.4 * 1e6 * 2*np.pi
    g = 420 * 2 * np.pi*1e6 

###

phase_vel=1.19*10**8
Z0 = 65.5
Zline = 50

Z_env = Zline
#print('eff_Cg:', eff_Cg)

omegas = np.linspace(5, 11, 250)* 2 * np.pi*1e9

## determine circuit parameters that exactly reproduce the frequencies of the coupled resonator system - without the notch

L_r, L_p, C_r, C_p, C_g, C_J, C_k = find_circuit_params_without_notch(freq_q, freq_r, freq_p, g, J, kappa_p, Cq, phase_vel, Z0, Z_env)

T1_without_notch = T1_limit_without_notch(omegas, L_r, L_p, C_r, C_p, C_g, C_J, C_k)
T1_without_notch = T1_without_notch*1e3

## determine circuit parameters that exactly reproduce the frequencies of the coupled resonator system - with the notch

L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k = find_circuit_params_with_notch(freq_q, freq_r, freq_p, freq_n, g, J, kappa_p, Cq, phase_vel, Z0, Z_env)

T1_with_notch = T1_limit_with_notch(omegas, L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k)
T1_with_notch = T1_with_notch*1e3

T1_with_notch_no_shunt = T1_limit_with_notch(omegas, L_r, L_p, C_r, C_p, C_g, L_J, C_J, C_k, with_shunt = False)
T1_with_notch_no_shunt = T1_with_notch_no_shunt*1e3

## determine circuit parameters that exactly reproduce the frequencies of the coupled resonator system - without the filter

# kappa_ro, _ =  hybridized_modes_kappas(freq_r, freq_p, kappa_p, J, chi, state = 'g')

# kappa_ro = kappa_ro #* 2

print('kappa_ro (MHz):', kappa_ro/ (2*np.pi*1e6))

T1_enhancement_val = T1_enhancement_factor(freq_q, freq_r, freq_p, freq_n, rounded = True)

print('T1_enhancement_val:', T1_enhancement_val)

T1_enhancement_from_filter_val = T1_enhancement_from_filter_factor(freq_q, freq_r, freq_p, J, kappa_p, kappa_ro)

print('T1_enhancement_from_filter_val:', T1_enhancement_from_filter_val)

T1_total_enhacement_val = T1_total_enhancement_factor(freq_q, freq_r, freq_p, freq_n, J, kappa_p, kappa_ro)

print('T1_total_enhacement_val (x10^3):', T1_total_enhacement_val/1e3)

L_r, C_r, C_g, C_k = find_circuit_params_without_filter(freq_q, freq_r, g, kappa_ro, Cq, phase_vel, Z0, Z_env)

T1_without_filter = T1_limit_without_filter(omegas, L_r, C_r, C_g, C_k)
T1_without_filter = T1_without_filter*1e3

T1_enhancement_vals = T1_enhancement_factor(omegas, freq_r, freq_p, freq_n, rounded = True)
T1_enhancement_from_filter_vals = T1_enhancement_from_filter_factor(omegas, freq_r, freq_p, J, kappa_p, kappa_ro)

plt.plot(omegas/(2*np.pi*1e9), T1_without_notch, label = 'without notch', color = 'blue')
plt.plot(omegas/(2*np.pi*1e9), T1_without_filter*T1_enhancement_from_filter_vals, label = 'without notch - analytic', color = 'blue', linestyle = '--')
plt.plot(omegas/(2*np.pi*1e9), T1_without_notch*T1_enhancement_vals, label = 'with notch - analytic', color = 'red', linestyle = '--')
plt.plot(omegas/(2*np.pi*1e9), T1_with_notch, label = 'with notch', color = 'red')
plt.plot(omegas/(2*np.pi*1e9), T1_without_filter, label = 'without filter', color = 'green')
plt.legend()
plt.xlabel('Frequency (GHz)')
plt.ylabel('T1 (ms)')
plt.title('T1 vs frequency')
plt.yscale('log')
plt.show()

# delta_qn = np.abs(omegas - freq_n)
# scaling = 1/4 * omegas **2 / (delta_qn**2) * (1-(freq_n/freq_r)**2)**2

# T1_without_notch_from_BW_eq = T1_with_notch_circ / scaling

# enhancement_at_qubit_freq = 1/4 * freq_q **2 / ((freq_q - freq_n)**2) * (1-(freq_n/freq_r)**2)**2
