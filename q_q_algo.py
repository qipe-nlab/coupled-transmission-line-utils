import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as opt

import q_q_formulas as q_q
import q_q_cap as cap
import common_formulas as cf
import q_r_cap

def function(inputs, target, f_q1, f_q2, res_len):
    angle = inputs[0]
    om_q1 = f_q1*2*np.pi

    Cq1 = 54e-15
    Cq2 = 45e-15

    Lq1 = 1/(Cq1*om_q1**2)
    Lq2 = Lq1
    print(Lq1)
    om_q2 = 1/np.sqrt(Lq2*Cq2)
    print("f2= ", om_q2/2/np.pi/1e9)
    Cc1 = cap.get_Cc(160e-6, angle)
    Cc2 = cap.get_Cc(140e-6, angle)

    Cr, Lr = cf.lumped_l_2_resonator_C_and_L(res_len)

    Z1 = q_q.lumped_model_Z21_no_ind(om_q1, Cq1, Cc1, Cr, Lr, Cc2, Cq2)
    Z2 = q_q.lumped_model_Z21_no_ind(om_q2, Cq1, Cc1, Cr, Lr, Cc2, Cq2)
    
    J = cf.lumped_elements_j_formula(om_q1, om_q2, Z1, Z2, Lq1, Lq2)
    
    return [J[0]-target[0]]


def print_message(sol, target, f_q1, f_q2, res_len):
    solution = sol.x
    result = function(solution, [0], f_q1, f_q2, res_len)
    
    print("\n \n")
    print("################################")
    print("     q_q: Solution found")
    print(" angle = ", int(solution[0]*10)/10, "[deg]")
    print("-------------------------------")
    print("         Results in")
    print(" J/2pi = ", int(result[0]*1e-3/2/np.pi)/1000, "[MHz]")
    print("################################")


def solve_for_q_q(target, x0, f_q1, fq_2, res_len, show=0):
    
    sol = opt.root(function, x0, args=(target, f_q1, f_q2, res_len))
    if show!=0:
        print_message(sol, target, f_q1, f_q2, res_len)
    return sol


f_q1 = 9e9
Lq1 = 5.79e-9

f_q2 = 10e9
Lq2 = 5.79e-9

res_len = 2255e-6


Z1s = np.array([0.59, 0.97, 1.35])*-1j
Z2s = np.array([0.56, 0.93, 1.29])*-1j

targets = np.array([10, 15, 20])
sim = np.zeros(3)
for i in range(3):
    sim[i] = cf.lumped_elements_j_formula(9e9*2*np.pi, 9.8e9*2*np.pi, Z1s[i], Z2s[i], Lq1, Lq2)/2/np.pi/1e6
    #print(sim[i])
unit = np.linspace(10,30,100)
plt.plot(unit, unit, linestyle="--", color="grey")
plt.plot(targets, sim, linestyle="none", marker="o", color="red")
plt.xlabel("targeted J [MHz]")
plt.ylabel("simulated J [MHz]")
plt.show()


"""

for i in [10,15,20]:
    x0 = [22.5]
    target=[i*1e6*2*np.pi]
    solve_for_q_q(target,x0,f_q1, f_q2, res_len, show=1)

"""
