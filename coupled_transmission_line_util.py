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
########################################

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
