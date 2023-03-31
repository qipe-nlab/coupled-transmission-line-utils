import numpy as np
import matplotlib.pyplot as plt
import sys

##########

### Basic functions

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

    val = (Lr/Cr)**0.5

    return val

def omega_r(Cr, Lr):

    val = 1/(Lr*Cr)**0.5

    return val

def transmission_line_voltage_current_out(Z0, electric_length, voltage_in, current_in):

    voltage_out = voltage_in*np.cos(electric_length)-1j*current_in*Z0*np.sin(electric_length)

    current_out = -1j*voltage_in/Z0*np.sin(electric_length)+current_in*np.cos(electric_length)

    return voltage_out, current_out

##########

### Weakly coupled lines model functions

def A_val(k_plus, k_minus, tau, gamma_Gf, omega):

    val = k_plus/(2*tau) * (1-np.exp(-1j*omega*2*tau)) - gamma_Gf * k_minus * np.exp(-1j*omega*2*tau) * 1j*omega 

    return val

def B_val(k_plus, k_minus, tau, gamma_Gf, omega):

    val = gamma_Gf * k_plus/(2*tau) * (1-np.exp(-1j*omega*2*tau)) * np.exp(-1j*omega*tau) - k_minus * 1j*omega * np.exp(-1j*omega*tau)

    return val

def reflection_ceofficient(ZL, Z0):

    val = (ZL - Z0)/(ZL + Z0)

    return val

def F_val(gamma1, gamma2, tau, omega):

    val = 1/(1-gamma1*gamma2*np.exp(-1j*omega*2*tau))

    return val

def Z_short_tl(Z0, phase_vel, L, omega):

    tau = L/phase_vel

    val = 1j*Z0*np.tan(tau * omega)

    return val 

def Z_open_tl(Z0, phase_vel, L, omega):

    tau = L/phase_vel

    val = -1j*Z0/np.tan(tau * omega)

    return val

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

def Z_transfer_total(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega):

    L_G = l_Gn + l_c + l_Gf

    #L_G = l_c + l_Gf

    #Z_trans_along_shorted_tr_val = Z_trans_along_shorted_tl(Z0, phase_vel, L_G, l_Gn, omega)

    Z_trans_along_shorted_tr_val = voltage_at_source_location(Z0, phase_vel, Cm_per_len, l_c, l_Gf, l_Gn, omega)

    voltage_transmission_coupled_lines_val = voltage_transmission_coupled_lines(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega)

    tau_Rn = l_Rn/phase_vel

    Rn_tl_voltage_scale_factor = 1/np.cos(tau_Rn*omega)

    val = Z_trans_along_shorted_tr_val * voltage_transmission_coupled_lines_val * Rn_tl_voltage_scale_factor

    return val

def lamda_by_4_omega(phase_vel, L):

    lambda_line = 4 * L

    val = 2 * np.pi * phase_vel / lambda_line

    return val

##########

Z0 = 50
phase_vel = 3*1e8/3
l_c = 30e-6 * 2
l_Gf = 400e-6 * 2
l_Gn = 900e-6 * 2
l_Rf = 900e-6 * 2
l_Rn = 500e-6 * 2
# Lm = 10000e-9
Cm_per_len = 50e-12
Lm_per_len = Z0**2 * Cm_per_len * 0.45

##########

### Lumped element model functions

def transmission_line_C_and_L(phase_vel, Z0):

    C_val = 1/(phase_vel*Z0)

    L_val = Z0 / phase_vel

    #print('C_val:', C_val)
    #print('L_val:', L_val)

    return C_val, L_val

def lumped_model_C_and_L(phase_vel, Z0, cpw__length):

    tl_C_val, tl_L_val = transmission_line_C_and_L(phase_vel, Z0)

    C_val = (tl_C_val * cpw__length / 2)*(0.9e-3+1e-3*(3e-15/1e-3+tl_C_val)/tl_C_val+1.7e-3)/cpw__length

    L_val = 8 * tl_L_val * cpw__length / np.pi**2 

    return C_val, L_val

def lumped_model_Cg_and_Lg(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm):

    omega_f = find_omega_f_analytic(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, 1e9*2*np.pi, 10e9*2*np.pi, 5e4*2*np.pi)

    Z0_f = find_filter_char_impedance(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f)
 
    Cg = 1/(omega_f*Z0_f)

    Lg =  Z0_f / omega_f

    return Cg, Lg

def lumped_model_transmission(L1, C1, L2, C2, Cg, Lg, omega):

    Z_res1 = Zres(C1, L1, omega)
    Z_res2 = Zres(C2, L2, omega)
    Z_res_coupler = Zres(Cg, Lg, omega)

    val = Z_res1 * Z_res2 / Z_res_coupler * 1/(1 + (Z_res1 + Z_res2)/Z_res_coupler)

    return val 

def lumped_model_resonator_coupling(L1, C1, L2, C2, Cg, Lg):

    omega_1 = omega_r(C1, L1)
    omega_2 = omega_r(C2, L2)

    Z21_omega_1 = lumped_model_transmission(1, C1, 1, C2, Cg, Lg, omega_1)
    Z21_omega_1 = lumped_model_transmission(1, C1, 1, C2, Cg, Lg, omega_2)

    val = -0.25*(omega_1**3*omega_2**3*C1*C2)**0.5*(np.imag(Z21_omega_1)/omega_1 + np.imag(Z21_omega_1)/omega_2)

    return val

##########

### find filter frequency and impedance

def find_omega_f_analytic(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, min_search, max_search, search_spacing):

    omegas = np.arange(min_search, max_search, search_spacing)

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

    results = defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas) - defining_eq2(Z0, phase_vel, Lm, Cm)

    # plt.plot(omegas/(1e9*2*np.pi), results)
    # plt.show()

    asign = np.sign(results)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)

    signchange[0] = 0

    # plt.plot(omegas/(1e9*2*np.pi), signchange)
    # plt.show()
    idxs = np.array(np.nonzero(signchange))[0]

    ### added debugger - check for continuity of eq1 around idx:
    print('idxs:', idxs)
    gaps = defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas[idxs + 1]) - defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas[idxs - 1])
    idx = idxs[(abs(gaps) < 1)]
    print('idxs:', idx)
    if idx.size == 0:
        raise ValueError('No valid solution to notch frequency equation for given input parameters in specified frequency range. Therefore cannot proceed to finding Lg and Cg.')

    
    Z_transfer_vals = np.abs(Z_transfer_total(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas[idx]))

    # print('Z_transfer_vals:', Z_transfer_vals)

    min_idx = np.argmin(Z_transfer_vals)

    idx = idx[min_idx]

    # print('idx:', idx)

    omega_f_rough  = omegas[idx]

    # print('omega_f_rough:', omega_f_rough/(1e9*2*np.pi))

    fine_search_spacing = 100 * 2*np.pi

    omegas = np.arange(omega_f_rough - search_spacing, omega_f_rough + search_spacing, fine_search_spacing)

    results = defining_eq1(Z0, phase_vel, Cm, l_c, l_Gf, l_Rf, omegas) - defining_eq2(Z0, phase_vel, Lm, Cm)

    asign = np.sign(results)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    signchange[0] = 0

    idx= np.nonzero(signchange)

    val  = omegas[idx][0]

    # plt.plot(omegas/(1e9*2*np.pi), signchange)
    # plt.show()
    # print (signchange)

    return val

def find_filter_char_impedance(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f):

    cpw__length1 = l_c + l_Gf + l_Gn
    cpw__length2 = l_c + l_Rf + l_Rn

    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1)
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2)

    Zres_1_omega_f = Zres(C1, L1, omega_f)
    Zres_2_omega_f = Zres(C2, L2, omega_f)

    delta_omega = 10 * 2*np.pi # 10 Hz

    Z_transfer_plus = Z_transfer_total(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f + delta_omega/2)
    Z_transfer_minus = Z_transfer_total(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omega_f - delta_omega/2)

    grad_Ztransfer = (Z_transfer_plus - Z_transfer_minus)/(delta_omega)

    val = np.real(1j*2*Zres_1_omega_f * Zres_2_omega_f / (omega_f * grad_Ztransfer))

    return val

##########

def find_notch_freq(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm):
    phase_vel = 3*10**8/2.5
    Z0 = 65
    return find_omega_f_analytic(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, 3*2*np.pi*10**9, 12*2*np.pi*10**9, 0.01*2*np.pi*10**9)/(2*np.pi*1e9)

def find_lambda_freq(cpw_length):
    phase_vel = 3*10**8/2.5
    return lamda_by_4_omega(phase_vel, cpw_length) / (2*np.pi*1e9)

def get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len):
    phase_vel = 3*10**8/2.5
    Z0 = 65
    cpw__length1 = l_c + l_Gf + l_Gn
    cpw__length2 = l_c + l_Rf + l_Rn
    
    C1, L1 = lumped_model_C_and_L(phase_vel, Z0, cpw__length1)
    C2, L2 = lumped_model_C_and_L(phase_vel, Z0, cpw__length2)
    Cg, Lg = lumped_model_Cg_and_Lg(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)

    Cg_basic = Cm_per_len *l_c
    Lg_basic = L1*L2/ (Lm_per_len * l_c)

    return C1, L1, C2, L2, Cg, Lg



# omega_vals = np.linspace(0.1, 15, 500) * 2*np.pi * 1e9

# Z_transfer_total_vals = Z_transfer_total(phase_vel, Z0, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len, omega_vals)



