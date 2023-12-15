import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
sys.path.append('./r_l')
sys.path.append('./r_r')
sys.path.append('./q_r')

import common_formulas as cf

import r_l_algo as r_l
import r_r_algo as r_r
import q_r_algo as q_r

import r_r_formulas as r_r_f

import cap_util as cap

def solve_all(f_r1, f_r2, f_notch, J, f_q, g, kappa, show=0):

    #
    # r_l
    #

    RL = 34
    r_l_x0 = [70, 2813e-6] # [Ck angle, len_r2]
    r_l_target = [kappa, f_r2]
    
    sol1 = r_l.solve_for_r_l(r_l_target, r_l_x0, RL, show=1)
    Ck_angle, len_r2 = sol1.x

    Ck = cap.get_Ck(Ck_angle)
    Cs = cf.calc_Cs(f_r2*2*np.pi, Ck, RL)
    Rs = cf.calc_Rs(f_r2*2*np.pi, Ck, RL)

    print(Rs)
    print(Cs)
    
    #
    # r_r
    #
    
    r_r_x0 = [1300e-6, 300e-6, 1300e-6, 15e-6] # [l_Gn, l_Gf, l_c, l_Rn, l_Rf, d]
    calibration_len=400e-6
    r_r_target = [cf.lambda_by_4_Ltot(f_r1), len_r2, f_notch, J]
    
    l_Gn, l_Gf, l_c, l_Rn, l_Rf, d = r_r.solve_for_r_r(r_r_target, r_r_x0, Cs)

    
    #
    # q_r
    #
    L_q = 5.5e-9
    q_l_x0 = [152e-6,30e-6] # [R_q, Cc_offset]
    q_l_target = [f_q, g]
    
    sol3 = q_r.solve_for_q_r(q_l_target, q_l_x0, f_r1, L_q, show=0)
    R_q, Cc_deformation = sol3.x
    #return 0

    if show==1:
        res_len_r2 = l_Rn+l_c+l_Rf+calibration_len
        
        # r-l calculations
        res_kappa, res_fr2 = r_l.function([Ck_angle, res_len_r2], [0,0], RL)


        # r-r calculations
        res_fr1 = cf.lambda_by_4_f(l_Gn+l_Gf+l_c)
        Cm_per_len = cap.get_Cm(d)
        Lm_per_len = cap.get_Lm(d)

        res_fn = r_r_f.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
        
        C1, L1, C2, L2, Cg, Lg = r_r_f.get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)
        res_j = r_r_f.lumped_model_get_j(C1, L1, C2+Cs, L2, Cg, Lg)

        # q-r calculations
        res_fq, res_g = q_r.function([R_q, Cc_deformation], [0,0], res_fr1, L_q)
        
        
        print("\n \n")
        print("################################")
        print("    Solution found")
        print("-------------------------------")
        print(" R_q = ", int(R_q*1e7)/10, "[um]")
        print(" Cc_deformation = ", int(Cc_deformation*1e7)/10, "[um]")
        print("-------------------------------")
        print(" Lgn = ", int(l_Gn*1e7)/10, "[um]")
        print(" Lgf = ", int(l_Gf*1e7)/10, "[um]")
        print(" Lc = ", int(l_c*1e7)/10, "[um]")
        print(" Lrn = ", int(l_Rn*1e7)/10, "[um]")
        print(" Lrf = ", int(l_Rf*1e7)/10, "[um]")
        print(" d = ", int(d*1e7)/10, "[um]")
        print("-------------------------------")
        print(" Ck_angle = ", int(Ck_angle*10)/10, " [deg]")
        print("-------------------------------")
        print("-------------------------------")
        print("         Results in")
        print(" f_q = ", int(res_fq*1e-6)/1000, "[GHz]")
        print(" g = ", int(res_g*1e-3/2/np.pi)/1000, "[MHz]")
        print("-------------------------------")
        print(" fr1 = ", int(res_fr1*1e-6)/1000, "[GHz]")
        print(" fr2 = ", int(res_fr2*1e-6)/1000, "[GHz]")
        print(" fn = ", int(res_fn*1e-6)/1000, "[GHz]")
        print(" J = ", int(res_j*1e-3/2/np.pi)/1000, "[MHz]")
        print("-------------------------------")
        print(" kappa/2pi = ", int(res_kappa/2/np.pi*1e-3)/1000, "[MHz]")
        print("################################")

        #print("To fix:")
        #r_r.get_notch_frequency_derivative(l_c, l_Gf, l_Gn, l_Rf, l_Rn, d, 50e-6)
        
    return [Ck_angle, l_Gn, l_Gf, l_c, l_Rn, l_Rf, d, R_q, Cc_deformation]

# Luka example parameters

# f_r1 = 6e9
# f_r2 = f_r1
# f_notch = 4e9
# J = 30e6 *2*np.pi

# f_q = f_notch
# g = 345e6 *2*np.pi

# kappa = 100e6 *2*np.pi

# new params

f_r1 = 10e9
f_r2 = f_r1
f_notch = 8e9
J = 5e6 *2*np.pi

f_q = f_notch
g = 100e6 *2*np.pi

kappa = 10e6 *2*np.pi

solve_all(f_r1, f_r2, f_notch, J, f_q, g, kappa, show=1)
