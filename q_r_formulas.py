
import numpy as np
import matplotlib.pyplot as plt
import sys

import common_formulas as cf


def lumped_model_Z21_no_ind(omega, C_q, L_q, C_r, L_r, C_g):
    Zq = 1/(1j*omega*C_q)
    Zg = 1/(1j*omega*C_g)
    Zr = 1/(1j*omega*C_r)

    Z_tot = Zq * Zr / Zg * 1/(1 + (Zq + Zr)/Zg)
    return Z_tot

def lumped_elements_get_coupling(omq, omr, Cq, Lq, cpw_len, Cg):
    Cr, Lr = cf.lumped_resonator_C_and_L(cpw_len)
    Z1 = lumped_model_Z21_no_ind(omq, Cq, Lq, Cr, Lr, Cg)
    Z2 = lumped_model_Z21_no_ind(omr, Cq, Lq, Cr, Lr, Cg)
    coupling = cf.lumped_elements_j_formula(omq, omr, Z1, Z2, Lq, Lr)
    return coupling

def g2r(g, om_q, om_r, c_q):
    return 2*g/(np.sqrt(om_q*om_r*c_q))

def r2g(r, om_q, om_r, c_q):
    return r*np.sqrt(om_q*om_r)/2*np.sqrt(c_q)


"""
#TESTING

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

res_len = cf.lambda_by_4_Ltot(om_r/(2*np.pi))

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
"""

r = 2.2e5
om_q = 6e9*2*np.pi
om_r = 8e9*2*np.pi
c_q = 6e-14

print(r2g(r, om_q, om_r, c_q)/1e6/2/np.pi)
