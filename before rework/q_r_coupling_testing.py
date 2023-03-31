import scipy as sp
import numpy as np
import matplotlib.pyplot as plt


from q_r_cap_util import get_Cq, get_Cg
#import common_formulas as cf


####################################################
### Taken from coupled_transmission_line_util.py ###
####################################################

# If needed in future: put into one shared library

def transmission_line_C_and_L(phase_vel, Z0):

    C_val = 1/(phase_vel*Z0)
    L_val = Z0 / phase_vel

    return C_val, L_val


def lumped_model_C_and_L(phase_vel, Z0, cpw__length):

    tl_C_val, tl_L_val = transmission_line_C_and_L(phase_vel, Z0)

    C_val = tl_C_val * cpw__length / 2
    L_val = 8 * tl_L_val * cpw__length / np.pi**2 

    return C_val, L_val

def lumped_elements_coupling_formula(om1, om2, Z1, Z2, L1, L2):
    return -0.25*np.sqrt(om1*om2/(L1*L2))*np.imag(Z1/om1+Z2/om2)

def q_r_Z21_no_ind(omega, C_q, L_q, C_r, L_r, C_g):
    Z1 = 1/(1j*omega*C_q)
    Z2 = 1/(1j*omega*C_g)
    Z3 = 1/(1j*omega*C_r)

    Z_tot = Z1 * Z3 / Z2 * 1/(1 + (Z1 + Z3)/Z2)
    return Z_tot

def lumped_elements_get_coupling(omq, omr, Cq, Lq, cpw_len, Cg):
    Cr, Lr = lumped_model_C_and_L(119919602, 65, cpw_len)
    Z1 = q_r_Z21_no_ind(omq, Cq, Lq, Cr, Lr, Cg)
    Z2 = q_r_Z21_no_ind(omr, Cq, Lq, Cr, Lr, Cg)
    coupling = lumped_elements_coupling_formula(omq, omr, Z1, Z2, Lq, Lr)
    return coupling

def lambda_by_4_omega(L, phase_vel=119919602):

    lambda_line = 4 * L
    omega = 2 * np.pi * phase_vel / lambda_line

    return omega


def lambda_by_4_Ltot(f, phase_vel=119919602):
    L = (phase_vel / f) / 4
    return L

################
### Formulas ###
################

def g2r(g, om_q, om_r, c_q):
    return 2*g/(np.sqrt(om_q*om_r*c_q))

def r2g(r, om_q, om_r, c_q):
    return r*np.sqrt(om_q*om_r)/2*np.sqrt(c_q)



om_r = 10.027e9*2*np.pi
r_sim = 3.01e5
R_q = 168e-6
offset = 37.8e-6

# Q1
#C_q = 43.62e-15
#C_q = 52.6e-15
#C_g = 7.78e-15

# Q2
#C_q = 34.3e-15
#C_g = 7.331e-15
#om_r = 10.640e9*2*np.pi

C_q = get_Cq(R_q)
C_g = get_Cg(offset)

res_len = lambda_by_4_Ltot(om_r/(2*np.pi))

f_q = 9e9
om_q = f_q*2*np.pi
L_q = 1/(om_q**2*C_q)


g = lumped_elements_get_coupling(om_q, om_r, C_q, L_q, res_len, C_g) # this assumes that J=g
r = g2r(g, om_q, om_r, C_q)

g_sim = r2g(r_sim, om_q, om_r, C_q)


print("qubit freq = ", 1/(np.sqrt(L_q*C_q))/2/np.pi/1e9)
print("resonator freq = ", om_r/1e9/2/np.pi)
print("resonator length = ", res_len*1e6)
print("g_calc/2pi = ", g/1e6/2/np.pi)
print("g_sim/2pi = ", g_sim/1e6/2/np.pi)
print("r_calc = ", r/1e5)
print("L_q = ", L_q)

