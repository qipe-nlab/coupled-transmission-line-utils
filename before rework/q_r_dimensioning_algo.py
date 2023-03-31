import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as opt

from q_r_coupling_testing import lumped_elements_get_coupling, lambda_by_4_Ltot
from q_r_cap_util import get_Cq, get_Cg

def function(inputs, target, f_r, L_q):
    R_q = inputs[0]
    coupler_offset = inputs[1]
    
    om_r = f_r*2*np.pi
    
    C_q = get_Cq(R_q)
    C_g = get_Cg(coupler_offset)

    res_len = lambda_by_4_Ltot(f_r)

    om_q = 1/(np.sqrt(L_q*C_q))
    f_q = om_q/2/np.pi
    
    g = lumped_elements_get_coupling(om_q, om_r, C_q, L_q, res_len, C_g)
    
    return [f_q-target[0], g-target[1]]


def print_message(sol, target, f_r, L_q):
    solution = sol.x
    result = function(solution, [0,0], f_r, L_q)
    
    print("\n \n")
    print("################################")
    print("       Solution found")
    print(" R_q = ", int(solution[0]*1e7)/10, "[um]")
    print(" coupler_offset = ", int(solution[1]*1e7)/10, "[um]")
    print("-------------------------------")
    print("         Results in")
    print(" f_q = ", int(result[0]*1e-6)/1000, "[GHz]")
    print(" g = ", int(result[1]*1e-3/2/np.pi)/1000, "[MHz]")
    print("################################")


def get_qubit_geometry(target, x0, f_r, L_q, show=0):
    
    sol = opt.root(function, x0, args=(target, f_r, L_q))
    if show!=0:
        print_message(sol, target, f_r, L_q)
    return sol


f_r = 10.027e9
L_q = 5.95e-9
    
x0 = [152e-6,30e-6]
target=[9e9, 320e6*2*np.pi] # target format: f_q, g
get_qubit_geometry(target,x0,f_r,L_q,show=1)
