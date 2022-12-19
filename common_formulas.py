import numpy as np
import matplotlib.pyplot as plt
import sys

"""
This file contains all the formulas that are used by multiple different parts of the dimensioning algorithm.

"""

# Global defaults
default_phase_vel = 119919602
default_Z0 = 65


#############################
# Lumped element impedances #
#############################

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


#################
# L/4 resonator #
#################

def lambda_by_4_omega(L, phase_vel=default_phase_vel):
    lambda_line = 4 * L
    val = 2 * np.pi * phase_vel / lambda_line
    return val

def lambda_by_4_f(L, phase_vel=default_phase_vel):
    lambda_line = 4 * L
    val = phase_vel / lambda_line
    return val

def lambda_by_4_Ltot(f, phase_vel=default_phase_vel):
    val = (phase_vel / f) / 4
    return val

def transmission_line_C_and_L(phase_vel=default_phase_vel, Z0=default_Z0):
    C_val = 1/(phase_vel*Z0)
    L_val = Z0 / phase_vel
    return C_val, L_val

def lumped_resonator_C_and_L(cpw__length, phase_vel=default_phase_vel, Z0=default_Z0):
    tl_C_val, tl_L_val = transmission_line_C_and_L(phase_vel, Z0)
    C_val = tl_C_val * cpw__length / 2
    L_val = 8 * tl_L_val * cpw__length / np.pi**2 

    return C_val, L_val

#################
# L/2 resonator #
#################

def lambda_by_2_omega(L, phase_vel=default_phase_vel):
    lambda_line = 2 * L
    val = 2 * np.pi * phase_vel / lambda_line
    return val

def lumped_l_2_resonator_C_and_L(cpw_length, phase_vel=default_phase_vel, Z0=default_Z0):
    omega = lambda_by_2_omega(cpw_length, phase_vel)
    C_val = np.pi/(2*omega*Z0)
    L_val = 1/(omega*C_val) 

    return C_val, L_val

######################
# Coupled resonators #
######################

def coupled_resonators_Z21(omega, Z1, Z2, Zc):
    return Z1 * Z2 / Zc * 1/(1 + (Z1 + Z2)/Zc)

def lumped_elements_j_formula(om1, om2, Z1, Z2, L1, L2):
    return -0.25*np.sqrt(om1*om2/(L1*L2))*np.imag(Z1/om1+Z2/om2)

####################
# Loaded resonator #
####################

def calc_Rs(om, Ck, RL):
    return (1+(om**2)*(Ck**2)*(RL**2))/((om**2)*(Ck**2)*RL)

def calc_Cs(om, Ck, RL):
    return Ck/(1+om**2*Ck**2*RL**2)
