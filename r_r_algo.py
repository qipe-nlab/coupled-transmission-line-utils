
import sys

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

import common_formulas as cf
import r_r_formulas as r_r
import cap_util as cap

def plot_transmission(l_Rf, l_Rn, l_Gf, l_Gn, l_c, d, Cs):

    omegas = np.linspace(3*2*np.pi*1e9, 10*2*np.pi*1e9, 10000)
    
    Cm_per_len = cap.get_Cm(d)
    Lm_per_len = cap.get_Lm(d)

    C1, L1, C2, L2, Cg, Lg = r_r.get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)
    Zs = r_r.lumped_model_transmission(C1, L1, C2, L2, Cg, Lg, omegas)
    Zs = np.abs(Zs)
    j_coupling = r_r.lumped_model_get_j(C1, L1, C2+Cs, L2, Cg, Lg)
    
    # notches
    notch_freq = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)

    
    print(f'notch_f = {notch_freq/2/np.pi*1e-9:.2f} GHz')
    print(f'J = {j_coupling/2/np.pi*1e-6:.2f} MHz')
    #   l/2
    cpw_length1 = l_c + l_Gf + l_Gn
    cpw_length2 = l_c + l_Rf + l_Rn
    cpw1_freq = cf.lambda_by_4_f(cpw_length1)
    cpw2_freq = cf.lambda_by_4_f(cpw_length2)
    

    plt.plot(omegas/2/np.pi/1e9, Zs)
    plt.vlines(x=notch_freq/2/np.pi/1e9, ymax=5, ymin=-5, color="red", label="notch_freq = " + str(int(notch_freq/1e6/2/np.pi)/1e3))
    plt.vlines(x=[cpw1_freq/1e9, cpw2_freq/1e9], ymax=5, ymin=-5, color="green")
    plt.xlabel("f [GHz]")
    plt.ylabel("Z21 [$ \Omega $]")
    plt.yscale('log')
    plt.legend()
    plt.show()

    

def solve_for_r_r(target, x0, Cs):
    """
    target format: [len1, len2, fn, J]
    x format: [Lgf, Lc, Lrf, d] (Lgn and Lrn will be determined based on the targeted lengths)
    """
    maxiter=1000
    len_stepsize=10e-6
    d_stepsize=0.2e-6
    reduced_len_stepsize=5e-6
    reduced_d_stepsize=0.1e-6
    
    f_tol=10e6
    J_tol=5e6
    
    state="notch"
    x=x0
    iter = 0
    success=0

    fn_log=[]
    J_log=[]
    state_log=[]

    
    while (iter<maxiter):
        l_Gn = target[0]-x[0]-x[1]
        l_Gf = x[0]
        l_c = x[1]
        l_Rn = target[1]-x[2]-x[1]
        l_Rf = x[2]
        d = x[3]

        
        Cm_per_len = cap.get_Cm(d)
        Lm_per_len = cap.get_Lm(d)

        notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
        
        C1, L1, C2, L2, Cg, Lg = r_r.get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)
        j_coupling = r_r.lumped_model_get_j(C1, L1, C2+Cs, L2, Cg, Lg)

        fn_log.append(notch_f/1e9)
        J_log.append(j_coupling/np.pi/2/1e6)
        state_log.append(state=="coupling")
        
        if iter%10==0:
            print("\r Current state:", state, " iteration:", iter, end="")
            sys.stdout.flush()

            # # enable for verbose debugging
            # print("\n\n")
            # print("iter = ", iter)
            
            # print(f'L_Gn = {l_Gn*1e6:.2f} um')
            # print(f'L_Gf = {l_Gf*1e6:.2f} um')
            # print(f'L_c = {l_c*1e6:.2f} um')
            # print(f'L_Rn = {l_Rn*1e6:.2f} um')
            # print(f'L_Rf = {l_Rf*1e6:.2f} um')
            # print(f'd = {d*1e6:.2f} um')
            
            # print(f'notch_f = {notch_f*1e-9:.2f} GHz')
            # print(f'J = {j_coupling/2/np.pi*1e-6:.2f} MHz')
            # print("state:", state)
            
        
        if state=="notch":
            if np.abs(target[2]-notch_f)<f_tol:
                state="coupling"
                #print("\n\n switch to coupling at ", iter)
            else:
                x[0]+=np.sign(notch_f-target[2])*len_stepsize
                x[2]+=np.sign(notch_f-target[2])*len_stepsize
            if np.abs(target[2]-notch_f)<2*f_tol and len_stepsize!=reduced_len_stepsize:
                #print("reduced len_stepsize")
                len_stepsize=reduced_len_stepsize

                
        if state=="coupling":
            if np.abs(target[3]-j_coupling)<J_tol:
                if np.abs(target[2]-notch_f)<f_tol:
                    state="finished"
                else:
                    #print("\n\n switching to notch at ", iter)
                    state="notch"
            else:
                if j_coupling>target[3] or d>=5e-6+d_stepsize:
                    x[3]+=np.sign(j_coupling-target[3])*d_stepsize
                else:
                    x[1]+=np.sign(target[3]-j_coupling)*len_stepsize
                    x[0]-=np.sign(target[3]-j_coupling)*len_stepsize/2
                    x[2]-=np.sign(target[3]-j_coupling)*len_stepsize/2
            if np.abs(target[3]-j_coupling)<2*J_tol and d_stepsize!=reduced_d_stepsize:
                #print("reduced d_stepsize")
                d_stepsize=reduced_d_stepsize

        if state=="finished":
            success = 1
            print("\n\n SUCCESS")
            break

        iter+=1
        
    plt.plot(range(len(state_log)), fn_log)
    plt.plot(range(len(state_log)), state_log)
    plt.show()
        
    print("----------")
    print("\n\n")

    if success==0:
        print("SEARCH FAILED")
    
    l_Gn = target[0]-x[0]-x[1]
    l_Gf = x[0]
    l_c = x[1]
    l_Rn = target[1]-x[2]-x[1]
    l_Rf = x[2]
    d = x[3]

    print(f'L_Gn = {l_Gn*1e6:.2f} um')
    print(f'L_Gf = {l_Gf*1e6:.2f} um')
    print(f'L_c = {l_c*1e6:.2f} um')
    print(f'L_Rn = {l_Rn*1e6:.2f} um')
    print(f'L_Rf = {l_Rf*1e6:.2f} um')
    print(f'd = {d*1e6:.2f} um')

    return [l_Gn, l_Gf, l_c, l_Rn, l_Rf, d]


"""
len_r1 = 5500e-6
len_r2 = 5300e-6
fn = 4e9
J = 20e6*2*np.pi
target = [len_r1, len_r2, fn, J]

Cs = 2.7e-14


x0 = [1300e-6, 300e-6, 1300e-6, 15e-6]
l_Gn, l_Gf, l_c, l_Rn, l_Rf, d = find_r_r(target, x0, Cs)

plot_transmission(l_Rf, l_Rn, l_Gf, l_Gn, l_c, d, Cs)
"""
