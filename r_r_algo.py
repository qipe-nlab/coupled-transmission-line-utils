import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

import common_formulas as cf
import r_r_formulas as r_r
import r_r_cap as cap


def function(inputs, target):
    calibration = 600e-6
    l_Gn = inputs[0]+calibration/2
    l_Gf = inputs[1]+calibration/2
    l_c = inputs[2]
    l_Rn = inputs[3]+calibration/2
    l_Rf = inputs[4]+calibration/2
    d = inputs[5]

    Cm_per_len = cap.get_Cm(d)
    Lm_per_len = cap.get_Lm(d)

    #
    # Resonator frequencies
    #
    f1 = cf.lambda_by_4_f(l_Gn+l_c+l_Gf)
    f2 = cf.lambda_by_4_f(l_Rn+l_c+l_Rf)

    #
    # Notch frequency
    #
    notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)

    #
    # J_coupling
    #
    j_coupling = r_r.lumped_model_get_j(f1*2*np.pi, f2*2*np.pi, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)

    #
    # Geometry constraints
    #
    constraint = 0
    
    # temporary constraint necessary for toy simulations
    # change when using different simulation model
    #if l_Gf/2>l_Gn-20e-6 or l_Rf/2>l_Rn-20e-6:
    #    constraint = 1e9

    if l_c>550e-6 or inputs[0]<800e-6 or inputs[3]<800e-6:
        constraint = 1e9
    
    return [f1-target[0], f2-target[1], notch_f-target[2], j_coupling-target[3], constraint, 0]


def print_message(sol, target, plot=0):
    solution = sol.x
    result = function(solution, [0,0,0,0])
    
    print("\n \n")
    print("################################")
    print("       Solution found")
    print(" Lgn = ", int(solution[0]*1e7)/10, "[um]")
    print(" Lgf = ", int(solution[1]*1e7)/10, "[um]")
    print(" Lc = ", int(solution[2]*1e7)/10, "[um]")
    print(" Lrn = ", int(solution[3]*1e7)/10, "[um]")
    print(" Lrf = ", int(solution[4]*1e7)/10, "[um]")
    print(" d = ", int(solution[5]*10)/10, "[um]")
    print("-------------------------------")
    print("         Results in")
    print(" fr1 = ", int(result[0]*1e-6)/1000, "[GHz]")
    print(" fr2 = ", int(result[1]*1e-6)/1000, "[GHz]")
    print(" fn = ", int(result[2]*1e-6)/1000, "[GHz]")
    print(" J = ", int(result[3]*1e-3)/1000, "[MHz]")
    print(" Constraints honored: ", result[4]==0)
    print("################################")

    if plot==1:
        omegas = np.linspace(10*np.pi*1e9, 24*np.pi*1e9, 10000)
        C1, L1, C2, L2, Cg, Lg = r_r.get_lumped_elements(solution[2], solution[1], solution[0], solution[4], solution[3], cap.get_Lm(solution[5]), cap.get_Cm(solution[5]),  phase_vel=119919602, Z0=64.6)
        Zs = r_r.lumped_model_Z21(omegas, C1, L1, C2, L2, Cg, Lg)
        Zs = np.abs(Zs)
        
        plt.plot(omegas/(2*np.pi*1e9), Zs, label="lumped element model", color="blue")
        plt.vlines([target[0]*1e-9,target[1]*1e-9], 1e-3, 1e5, color="orange", label="resonators")
        plt.vlines([target[2]*1e-9], 1e-3, 1e5, color="orange", label="notch")
        plt.yscale("log")
        plt.show()
    return 0

def get_readout_filter_geometry(target, x0, show=0, plot=0):
    sol = opt.root(function, x0, args=(target))
    if show!=0:
        print_message(sol, target, plot=plot)
    return sol




x0 = [800e-6, 1600e-6, 370e-6, 827e-6, 1100e-6, 5]
target=[9e9, 10e9, 8e9, 200e6] # target format: Fr1, Fr2, Fn, J
get_readout_filter_geometry(target,x0,1,1)


