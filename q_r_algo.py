import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as opt

import q_r_formulas as q_r
import q_r_cap as cap
import common_formulas as cf

def function(inputs, target, f_r, L_q):
    R_q = inputs[0]
    coupler_offset = inputs[1]
    #print("R_q = ", R_q*1e6)
    om_r = f_r*2*np.pi
    #print(R_q*1e6)
    C_q = cap.get_Cq(R_q)
    C_g = cap.get_Cg(coupler_offset)
    
    res_len = cf.lambda_by_4_Ltot(f_r)

    om_q = 1/(np.sqrt(L_q*C_q))
    f_q = om_q/2/np.pi
    #print("f_q = ", f_q/1e9)
    #print(f_q/1e9)
    g = q_r.lumped_elements_get_coupling(om_q, om_r, C_q, L_q, res_len, C_g)
    return [f_q-target[0], g-target[1]]


def print_message(sol, target, f_r, L_q):
    solution = sol.x
    result = function(solution, [0,0], f_r, L_q)
    
    print("\n \n")
    print("################################")
    print("     q_r: Solution found")
    print(" R_q = ", int(solution[0]*1e7)/10, "[um]")
    print(" coupler_offset = ", int(solution[1]*1e7)/10, "[um]")
    print("-------------------------------")
    print("         Results in")
    print(" f_q = ", int(result[0]*1e-6)/1000, "[GHz]")
    print(" g = ", int(result[1]*1e-3/2/np.pi)/1000, "[MHz]")
    print("################################")


def solve_for_q_r(target, x0, f_r, L_q, show=0):
    
    sol = opt.root(function, x0, args=(target, f_r, L_q))
    if show!=0:
        print_message(sol, target, f_r, L_q)
    return sol

"""
f_r = 9e9
L_q = 5.95e-9
    
x0 = [152e-6,30e-6]
target=[9e9, 260e6*2*np.pi] # target format: f_q, g
solve_for_q_r(target,x0,f_r,L_q,show=1)
"""
