import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

import cap_util as cap
import r_l_formulas as r_l
import common_formulas as cf



def function(inputs, target, RL):
    angle = inputs[0]
    len_r2 = inputs[1]
    
    Ck = cap.get_Ck(angle)
    om = cf.lambda_by_4_omega(len_r2)
    
    C, L = cf.lumped_resonator_C_and_L(len_r2)

    Rs = cf.calc_Rs(om, Ck, RL)
    Cs = cf.calc_Cs(om, Ck, RL)

    oms = cf.omega_r(L, C+Cs)
    Q_ext = oms*Rs*C

    fs = oms/2/np.pi
    
    kappa = oms/Q_ext
    return [kappa-target[0], fs-target[1]]

def print_message(sol, target, RL):
    solution = sol.x
    result = function(solution, [0,0], RL)
    
    print("\n \n")
    print("################################")
    print("    r_l: Solution found")
    print(" angle = ", int(solution[0]*10)/10, " [deg]")
    print(" len_r2 = ", int(solution[1]*1e7)/10, " [um]")
    print("-------------------------------")
    print("         Results in")
    print(" kappa/2pi = ", int(result[0]/2/np.pi*1e-3)/1000, "[MHz]")
    print(" f_r2 = ", int(result[1]*1e-6)/1000, "[GHz]")
    print("################################")


def solve_for_r_l(target, x0, RL, show=0):
    
    sol = opt.root(function, x0, args=(target, RL))
    if show!=0:
        print_message(sol, target, RL)
    return sol

"""
r = np.linspace(10, 10000, 1000)

RL = 34

angle_0 = 35
len_res = 2813e-6 #This would be the unloaded resonator length
#f_res = cf.lambda_by_4_f(len_res)
#print(f_res/1e9)

target = [40e6*2*np.pi, 10e9]

solve_for_r_l(target, [angle_0, len_res], RL, show=1)
"""
#print(function(angle_0, 0, f_res, RL)/2/np.pi/1e6)



#sim_k = np.array([9.6e6, 1.86e7, 2.9e7, 4.42e7, 6e7])
#angles = [10, 22.5, 35, 47.5, 60]
#calc_k = function(get_Ck(angles), 0, f_res, RL)/2/np.pi
"""
plt.plot(angles, sim_k/1e6, linestyle="none", marker="o")
plt.plot(angles[0], sim_k[0]/1e6, linestyle="none", marker="x", color="red", markersize=20)
plt.plot(angles, function(angles, 0, f_res, RL)/2/np.pi/1e6, linestyle="none", marker=".")
plt.xlabel("angle [deg]")
plt.ylabel("kappa/2pi [MHz]")
plt.show()
"""
