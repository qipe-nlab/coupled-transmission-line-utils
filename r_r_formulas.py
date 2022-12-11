import numpy as np
import matplotlib.pyplot as plt
import sys

import common_formulas as cf


########################################
### Distributed Element Calculations ###
### 2 coupled transmission lines, WC ###
########################################

### For 2 coupled transmission lines under the assumption of weak coupling
### See: Solution of the Transmission-Line Equations Under the Weak-Coupling Assumption

def transmission_line_voltage_current_out(electric_length, voltage_in, current_in, Z0=cf.default_Z0):
    voltage_out = voltage_in*np.cos(electric_length)-1j*current_in*Z0*np.sin(electric_length)
    current_out = -1j*voltage_in/Z0*np.sin(electric_length)+current_in*np.cos(electric_length)
    return voltage_out, current_out


def A_val(k_plus, k_minus, tau, gamma_Gf, omega):
    return k_plus/(2*tau) * (1-np.exp(-1j*omega*2*tau)) - gamma_Gf * k_minus * np.exp(-1j*omega*2*tau) * 1j*omega


def B_val(k_plus, k_minus, tau, gamma_Gf, omega):
    return gamma_Gf * k_plus/(2*tau) * (1-np.exp(-1j*omega*2*tau)) * np.exp(-1j*omega*tau) - k_minus * 1j*omega * np.exp(-1j*omega*tau)


def reflection_ceofficient(ZL, Z0=cf.default_Z0):
    return (ZL - Z0)/(ZL + Z0)


def F_val(gamma1, gamma2, tau, omega):
    return 1/(1-gamma1*gamma2*np.exp(-1j*omega*2*tau))


def Z_short_tl(L, omega, phase_vel=cf.default_phase_vel, Z0=cf.default_Z0):
    tau = L/phase_vel
    return 1j*Z0*np.tan(tau * omega)


def Z_open_tl(L, omega, phase_vel=cf.default_phase_vel, Z0=cf.default_Z0):
    tau = L/phase_vel
    return -1j*Z0/np.tan(tau * omega)


def voltage_transmission_coupled_lines(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega, phase_vel=cf.default_phase_vel, Z0=cf.default_Z0):

    Lm = Lm_per_len * l_c
    Cm = Cm_per_len * l_c

    #### testing using exact expressions:
    C_tl, L_tl = cf.transmission_line_C_and_L(phase_vel, Z0)
    Z0_c = cf.Zchar(C_tl + Cm_per_len, L_tl)
    phase_vel_c = cf.omega_r(C_tl + Cm_per_len, L_tl)

    K_plus = Lm/Z0_c + Cm*Z0_c

    K_minus = Lm/Z0_c - Cm*Z0_c

    tau_c = l_c / phase_vel_c

    Z_Rn = Z_open_tl(l_Rn, omega, phase_vel, Z0)
    
    Z_Rf = Z_short_tl(l_Rf, omega, phase_vel, Z0)
    
    Z_Gf = Z_short_tl(l_Gf, omega, phase_vel, Z0)

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


def Z_trans_along_shorted_tl(L, z, omega, phase_vel=cf.default_phase_vel, Z0=cf.default_Z0):
    tau = L/phase_vel
    val = 1j * Z0 * (np.tan(tau*omega)*np.cos(tau*z/L * omega) - np.sin(tau*z/L * omega)) 
    return val


def voltage_at_source_location(Cm_per_len, l_c, l_Gf, l_Gn, omega, phase_vel=cf.default_phase_vel, Z0=cf.default_Z0):
    Cl, Ll = cf.transmission_line_C_and_L(phase_vel, Z0)
    Z0_c = cf.Zchar(Cl + Cm_per_len, Ll)
    phase_vel_c = cf.omega_r(Cl + Cm_per_len, Ll)
    tau_f = l_Gf / phase_vel
    Zinput_f = cf.Zinput(Z0, 0, omega * tau_f)
    tau_c = l_c/phase_vel_c
    Zinput_c = cf.Zinput(Z0_c, Zinput_f, omega * tau_c)
    tau_n = l_Gn / phase_vel
    Zinput_n = cf.Zinput(Z0, Zinput_c, omega * tau_n)
    I_in = 1
    v_in = Zinput_n * I_in
    v_out, I_out = transmission_line_voltage_current_out(tau_n, v_in, I_in, Z0)
    val = v_out
    return val


def Z_transfer_total(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega, phase_vel=cf.default_phase_vel, Z0=cf.default_Z0):
    L_G = l_Gn + l_c + l_Gf
    #L_G = l_c + l_Gf
    #Z_trans_along_shorted_tr_val = Z_trans_along_shorted_tl(L_G, l_Gn, omega, phase_vel, Z0)
    Z_trans_along_shorted_tr_val = voltage_at_source_location(Cm_per_len, l_c, l_Gf, l_Gn, omega, phase_vel, Z0)
    voltage_transmission_coupled_lines_val = voltage_transmission_coupled_lines(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega, phase_vel, Z0)
    tau_Rn = l_Rn/phase_vel
    Rn_tl_voltage_scale_factor = 1/np.cos(tau_Rn*omega)
    val = Z_trans_along_shorted_tr_val * voltage_transmission_coupled_lines_val * Rn_tl_voltage_scale_factor
    return val


#1e9*2*np.pi, 10e9*2*np.pi, 5e4*2*np.pi
#min_search=3*2*np.pi*10**9, max_search=12*2*np.pi*10**9, search_spacing=(0.01*2*np.pi*10**9)
def find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=cf.default_phase_vel, Z0=cf.default_Z0, min_search=1e9*2*np.pi, max_search=10e9*2*np.pi, search_spacing=5e4*2*np.pi):

    # Define the two sides of the equation
    # We will find the zero crossing, to get the solution to the equation
    # Note: omega is the variable in which we want to find the crossing
    def defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omega):
        Cl, Ll= cf.transmission_line_C_and_L(phase_vel, Z0)
        phase_vel_c = cf.omega_r(Cl + Cm, Ll)
        tau_c = l_c / phase_vel_c
        tau_G_dash = (l_Gf / phase_vel + l_c / phase_vel_c / 2) 
        tau_R_dash = (l_Rf / phase_vel + l_c / phase_vel_c / 2) 
        val = omega * tau_c / (np.sin(omega * tau_c)) * (np.cos(2*omega*tau_G_dash) + np.cos(2*omega*tau_R_dash)) /(1 + np.cos(2*omega*(tau_G_dash + tau_R_dash)))
        return val

    def defining_eq2(Z0, phase_vel, Lm, Cm):
        Cl, Ll= cf.transmission_line_C_and_L(phase_vel, Z0)
        Zm = cf.Zchar(Cm, Lm)
        Z0_c = cf.Zchar(Cl + Cm, Ll)
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
    #print("idxs -> ", idxs)

    #plt.plot(omegas,  defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas) - defining_eq2(Z0, phase_vel, Lm, Cm), marker=".")
    #plt.plot(omegas[idxs], defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas[idxs]), color="red", marker="x")
    #plt.show()
    #quit()
    
    # added debugger - check for continuity of eq1 around idx:
    gaps = defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas[idxs + 1]) - defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas[idxs - 1])
    idx = idxs[(abs(gaps) < 1)]
    if idx.size == 0:
        #print('idxs:', idxs)
        raise ValueError('No valid solution to notch frequency equation for given input parameters in specified frequency range. Therefore cannot proceed to finding Lg and Cg.')

    
    Z_transfer_vals = np.abs(Z_transfer_total(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas[idx]))


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


def find_notch_filter_char_impedance(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f, phase_vel=cf.default_phase_vel, Z0=cf.default_Z0):

    cpw__length1 = l_c + l_Gf + l_Gn
    cpw__length2 = l_c + l_Rf + l_Rn

    C1, L1 = cf.lumped_resonator_C_and_L(cpw__length1, phase_vel, Z0)
    C2, L2 = cf.lumped_resonator_C_and_L(cpw__length2, phase_vel, Z0)

    Zres_1_omega_f = cf.Zres(C1, L1, omega_f)
    Zres_2_omega_f = cf.Zres(C2, L2, omega_f)

    delta_omega = 10 * 2*np.pi # 10 Hz

    Z_transfer_plus = Z_transfer_total(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f + delta_omega/2, phase_vel=phase_vel, Z0=Z0)
    Z_transfer_minus = Z_transfer_total(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f - delta_omega/2, phase_vel=phase_vel, Z0=Z0)

    grad_Ztransfer = (Z_transfer_plus - Z_transfer_minus)/(delta_omega)

    val = np.real(1j*2*Zres_1_omega_f * Zres_2_omega_f / (omega_f * grad_Ztransfer))
    
    return val


############################
### Lumped Element Model ###
############################

def lumped_model_Cg_and_Lg(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel=cf.default_phase_vel, Z0=cf.default_Z0):

    omega_f = find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0, 1e9*2*np.pi, 10e9*2*np.pi, 5e4*2*np.pi)
    Z0_f = find_notch_filter_char_impedance(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega_f, phase_vel=cf.default_phase_vel, Z0=Z0)

    Cg = 1/(omega_f*Z0_f)
    Lg =  Z0_f / omega_f

    return Cg, Lg


def get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel=cf.default_phase_vel, Z0=cf.default_Z0):
        
    cpw__length1 = l_c + l_Gf + l_Gn
    cpw__length2 = l_c + l_Rf + l_Rn
    
    C1, L1 = cf.lumped_resonator_C_and_L(cpw__length1, phase_vel, Z0)
    C2, L2 = cf.lumped_resonator_C_and_L(cpw__length2, phase_vel, Z0)
    Cg, Lg = lumped_model_Cg_and_Lg(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, phase_vel, Z0)
    
    #Cg_basic = Cm_per_len *l_c
    #Lg_basic = L1*L2 / (Lm_per_len * l_c)

    return C1, L1, C2, L2, Cg, Lg


def lumped_model_Z21(omega, C1, L1, C2, L2, Cg, Lg):
    Z1 = 1/(1j*omega*C1 + 1/(1j*omega*L1))
    Zg = 1/(1j*omega*Cg + 1/(1j*omega*Lg))
    Z2 = 1/(1j*omega*C2 + 1/(1j*omega*L2))

    Z_tot = cf.coupled_resonators_Z21(omega, Z1, Z2, Zg)
    return Z_tot


def lumped_model_Z21_no_ind(omega, C1, L1, C2, L2, Cg, Lg):
    Z1 = 1/(1j*omega*C1)
    Zg = 1/(1j*omega*Cg + 1/(1j*omega*Lg))
    Z2 = 1/(1j*omega*C2)

    Z_tot = cf.coupled_resonators_Z21(omega, Z1, Z2, Zg)
    return Z_tot


def lumped_model_get_j(om1, om2, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len):
    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)
    Z1 = lumped_model_Z21_no_ind(om1, C1, L1, C2, L2, Cg, Lg)
    Z2 = lumped_model_Z21_no_ind(om2, C1, L1, C2, L2, Cg, Lg)
    j = cf.lumped_elements_j_formula(om1, om2, Z1, Z2, L1, L2)
    return j
