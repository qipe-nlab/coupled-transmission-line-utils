import numpy as np
import matplotlib.pyplot as plt
import sys

import common_formulas as cf

def lumped_model_Z21_no_ind(omega, Cq1, Cc1, Cr, Lr, Cc2, Cq2):

    Zq1 = cf.Zcap(Cq1, omega)
    Zc1 = cf.Zcap(Cc1, omega)
    Zr = cf.Zres(Cr, Lr, omega)
    Zc2 = cf.Zcap(Cc2, omega)
    Zq2 = cf.Zcap(Cq2, omega)

    # All the following currents are scaled by V2
    Iq2 = 1/Zq2
    Vr = 1 + Iq2*Zc2
    Ir = Vr/Zr
    V1 = Vr + Zc1 * (Iq2+Ir)
    Iq1 = V1/Zq1

    # Z21 = V2/I1
    # See documentation for the derivation of this formula
    Z21 = 1/(Iq2 + Ir + Iq1)

    return Z21
"""
Cq1 = 53.8e-15
Cc1 = 5.08e-15
Cq2 = 40.6e-15
Cc2 = 4.6e-15

Lq = 5.95e-9

f1 = cf.omega_r(Cq1, Lq)/2/np.pi
f2 = cf.omega_r(Cq2, Lq)/2/np.pi

omega1 = f1*2*np.pi
omega2 = f2*2*np.pi


res_len = 2255e-6
Cr, Lr = cf.lumped_l_2_resonator_C_and_L(res_len)

print(np.abs(lumped_model_Z21_no_ind(omega1, Cq1, Cc1, Cr, Lr, Cc2, Cq2)))

Z1 = lumped_model_Z21_no_ind(omega1, Cq1, Cc1, Cr, Lr, Cc2, Cq2)
Z2 = lumped_model_Z21_no_ind(omega2, Cq1, Cc1, Cr, Lr, Cc2, Cq2)


J = cf.lumped_elements_j_formula(omega1, omega2, Z1, Z2, Lq, Lq)
    
print(J/2/np.pi/1e6)
"""
