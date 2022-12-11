import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

import r_l_cap as cap
import r_l_formulas as r_l
import common_formulas as cf



def function(angle, target, f_res, RL):
    Ck = cap.get_Ck(angle)
    om = 2*np.pi*f_res
    cpw_len = cf.lambda_by_4_Ltot(f_res)
    C, L = cf.lumped_resonator_C_and_L(cpw_len)

    Rs = r_l.calc_Rs(om, Ck, RL)
    Cs = r_l.calc_Cs(om, Ck, RL)

    oms = cf.omega_r(L, C+Ck)
    Q_ext = oms*Rs*C

    kappa = oms/Q_ext
    return kappa-target

def print_message(sol, target, f_r, RL):
    solution = sol.x
    result = function(solution, 0, f_r, RL)
    
    print("\n \n")
    print("################################")
    print("       Solution found")
    print(" angle = ", int(solution[0]*10)/10, " deg")
    print("-------------------------------")
    print("         Results in")
    print(" kappa/2pi = ", int(result/2/np.pi*1e-3)/1000, "[MHz]")
    print("################################")


def solve_for_angle(target, x0, f_r, RL, show=0):
    
    sol = opt.root(function, x0, args=(target, f_r, RL))
    if show!=0:
        print_message(sol, target, f_r, RL)
    return sol


r = np.linspace(10, 10000, 1000)

RL = 34
len_res = 2813e-6
f_res = cf.lambda_by_4_f(len_res)
print(f_res/1e9)
#f_res = 10.8e9

angle_0 = 35
target = 15e6*2*np.pi

solve_for_angle(target, angle_0, f_res, RL, show=1)

#print(function(angle_0, 0, f_res, RL)/2/np.pi/1e6)



sim_k = np.array([9.6e6, 1.86e7, 2.9e7, 4.42e7, 6e7])
angles = [10, 22.5, 35, 47.5, 60]
#calc_k = function(get_Ck(angles), 0, f_res, RL)/2/np.pi
"""
plt.plot(angles, sim_k/1e6, linestyle="none", marker="o")
plt.plot(angles[0], sim_k[0]/1e6, linestyle="none", marker="x", color="red", markersize=20)
plt.plot(angles, function(angles, 0, f_res, RL)/2/np.pi/1e6, linestyle="none", marker=".")
plt.xlabel("angle [deg]")
plt.ylabel("kappa/2pi [MHz]")
plt.show()
"""
