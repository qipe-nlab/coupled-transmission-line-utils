import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import scipy.optimize as opt
from r_v_cap_util import get_Ck


def transmission_line_C_and_L(phase_vel, Z0):

    C_val = 1/(phase_vel*Z0)
    L_val = Z0 / phase_vel

    return C_val, L_val

def lumped_model_C_and_L(phase_vel, Z0, cpw__length):

    tl_C_val, tl_L_val = transmission_line_C_and_L(phase_vel, Z0)

    C_val = tl_C_val * cpw__length / 2
    L_val = 8 * tl_L_val * cpw__length / np.pi**2 

    return C_val, L_val

def lambda_by_4_Ltot(f, phase_vel=119919602):
    L = (phase_vel / f) / 4
    return L

def lambda_by_4_f(L, phase_vel=119919602):
    lambda_line = 4 * L
    val = phase_vel / lambda_line
    return val


def function(angle, target, f_res, RL):
    Ck = get_Ck(angle)
    om = 2*np.pi*f_res
    cpw_len = lambda_by_4_Ltot(f_res)
    C, L = lumped_model_C_and_L(119919602, 65, cpw_len)

    Rs = (1+(om**2)*(Ck**2)*(RL**2))/((om**2)*(Ck**2)*RL)
    Cs = Ck/(1+om**2*Ck**2*RL**2)

    oms = 1/np.sqrt(L*(C+Cs))
    Q_ext = om*Rs*C

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
f_res = lambda_by_4_f(len_res)
print(f_res/1e9)
#f_res = 10.8e9

angle_0 = 35
target = 15e6*2*np.pi

#solve_for_angle(target, angle_0, f_res, RL, show=1)

#print(function(angle_0, 0, f_res, RL)/2/np.pi/1e6)



sim_k = np.array([9.6e6, 1.86e7, 2.9e7, 4.42e7, 6e7])
angles = [10, 22.5, 35, 47.5, 60]
#calc_k = function(get_Ck(angles), 0, f_res, RL)/2/np.pi


plt.plot(angles, sim_k/1e6, linestyle="none", marker="o", label="simulation")
plt.plot(angles, function(angles, 0, f_res, RL)/2/np.pi/1e6, linestyle="none", marker=".", label="formula")
plt.plot(angles[0], sim_k[0]/1e6, linestyle="none", marker="x", color="red", markersize=20, label="calibration point")
plt.legend()
plt.xlabel("angle [deg]")
plt.ylabel("kappa/2pi [MHz]")
plt.show()

