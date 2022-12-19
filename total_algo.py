import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


import common_formulas as cf
import r_l_algo as r_l
import r_r_algo as r_r
import q_r_algo as q_r

import r_l_cap as r_l_cap

def solve_all(f_r1, f_r2, f_notch, J, f_q, g, kappa):

    #
    # r_l
    #
    RL = 34
    r_l_x0 = [35, 2813e-6] # [Ck angle, len_r2]  
    r_l_target = [kappa, f_r2]
    
    sol1 = r_l.solve_for_r_l(r_l_target, r_l_x0, RL, show=1)
    Ck_angle, len_r2 = sol1.x

    Ck = r_l_cap.get_Ck(Ck_angle)
    Cs = cf.calc_Cs(f_r2*2*np.pi, Ck, RL)
    
    #
    # r_r
    #
    
    r_r_x0 = [1000e-6, 1000e-6, 1000e-6, 1000e-6, 1000e-6, 5e-6] # [l_Gn, l_Gf, l_c, l_Rn, l_Rf, d]
    r_r_x0 = [1400e-6, 1235e-6, 500e-6, 1400e-6, 1100e-6, 4e-6] # [l_Gn, l_Gf, l_c, l_Rn, l_Rf, d]
    r_r_target = [f_r1, len_r2, f_notch, J]
    
    sol2 = r_r.solve_for_r_r(r_r_target, r_r_x0, Cs, show=1, plot=0)
    l_Gn, l_Gf, l_c, l_Rn, l_Rf, d = sol2.x
    
    #
    # q_l
    #
    L_q = 15e-9
    q_l_x0 = [152e-6,30e-6] # [R_q, Cc_offset]
    q_l_target = [f_q, g]
    
    sol3 = q_r.solve_for_q_r(q_l_target, q_l_x0, f_r1, L_q, show=1)
    R_q, Cc_offset = sol3.x
    #return 0
    return [Ck_angle, l_Gn, l_Gf, l_c, l_Rn, l_Rf, d, R_q, Cc_offset]


f_r1 = 8e9
f_r2 = 8e9
f_notch = 6e9
J = 50e6 *2*np.pi


f_q = 6e9
g = 200e6 *2*np.pi

kappa = 10e6 *2*np.pi


solve_all(f_r1, f_r2, f_notch, J, f_q, g, kappa)
