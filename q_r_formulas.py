
import numpy as np
import matplotlib.pyplot as plt
import sys

import common_formulas as cf


def lumped_model_Z21_no_ind(omega, C_q, L_q, C_r, L_r, C_g):
    Zq = 1/(1j*omega*C_q)
    Zg = 1/(1j*omega*C_g)
    #Zr = 1/(1j*omega*C_r)
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

