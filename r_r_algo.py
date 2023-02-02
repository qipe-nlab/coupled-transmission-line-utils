import sys

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

import common_formulas as cf
import r_r_formulas as r_r
import cap_util as cap


def function(inputs, target, Cs, calibration_len, weighted_error=0):
    #print(inputs)
    l_Gn = inputs[0]+calibration_len/2
    l_Gf = inputs[1]+calibration_len/2
    l_c = inputs[2]
    l_Rn = inputs[3]+calibration_len/2
    l_Rf = inputs[4]+calibration_len/2
    d = inputs[5]

    Cm_per_len = cap.get_Cm(d)
    Lm_per_len = cap.get_Lm(d)

    #
    # Resonator frequencies
    #
    f1 = cf.lambda_by_4_f(l_Gn+l_c+l_Gf)
    #f2 = cf.lambda_by_4_f(l_Rn+l_c+l_Rf)
    len_r2 = l_Rn+l_c+l_Rf

    #
    # Notch frequency
    #
    notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)

    #
    # J_coupling
    #
    C1, L1, C2, L2, Cg, Lg = r_r.get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)
    om2 = cf.omega_r(C2+Cs, L2)
    #print(inputs)
    #print(C1, L1, C2, L2, Cg, Lg)
    j_coupling = r_r.lumped_model_get_j(C1, L1, C2+Cs, L2, Cg, Lg)
    omegas = np.linspace(10*np.pi*1e9, 24*np.pi*1e9, 10000)
    Zs = r_r.lumped_model_Z21(omegas, C1, L1, C2+Cs, L2, Cg, Lg)
    Zs = np.abs(Zs)
        
    #plt.plot(omegas/(2*np.pi*1e9), Zs, label="lumped element model", color="blue")
    #plt.yscale("log")
    #plt.show()

    
    #
    # Geometry constraints
    #
    constraint = 0
    
    # temporary constraint necessary for toy simulations
    # change when using different simulation model
    #if l_Gf/2>l_Gn-20e-6 or l_Rf/2>l_Rn-20e-6:
    #    constraint = 1e9

    #if l_c>850e-6 or inputs[0]<800e-6 or inputs[3]<800e-6:
    #    constraint = 1e12

    if inputs[0]<800e-6 or inputs[1]<10e-6 or inputs[2]<10e-6 or inputs[3]<800e-6 or inputs[4]<10e-6:
        constraint = 1e12

    #print("J = ", j_coupling/2/np.pi/1e6)
    #print("notch_f = ", notch_f/1e9)
    #print("f1", f1/1e9)
    #print("f2", cf.omega_r(C2+Cs, L2)/2/np.pi/1e9)
    #print("-------")
    #print([(f1-target[0])/1e6, (len_r2-target[1])*1e6, (notch_f-target[2])/1e6, (j_coupling-target[3])/1e6, constraint, 0])
    weights=[1,1,1,1,1]
    if weighted_error: return [(f1-target[0])*weights[0], (len_r2-target[1])*weights[1], (notch_f-target[2])*weights[2], (j_coupling-target[3])*weights[3], (constraint)*weights[4], 0]
    
    return [f1-target[0], len_r2-target[1], notch_f-target[2], j_coupling-target[3], constraint, 0]



def function2(inputs, target, Cs, calibration_len, weighted_error=1):
    l_Gn = inputs[0] +calibration_len/2
    l_Gf = inputs[1] +calibration_len/2
    l_c = inputs[2]
    l_Rn = inputs[3] +calibration_len/2
    l_Rf = inputs[4] +calibration_len/2
    d = inputs[5]
        
    Cm_per_len = cap.get_Cm(d)
    Lm_per_len = cap.get_Lm(d)

    #
    # Resonator frequencies
    #
    f1 = cf.lambda_by_4_f(l_Gn+l_c+l_Gf)
    #f2 = cf.lambda_by_4_f(l_Rn+l_c+l_Rf)
    len_r2 = l_Rn+l_c+l_Rf

    #
    # Notch frequency
    #
    notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)

    #
    # J_coupling
    #
    C1, L1, C2, L2, Cg, Lg = r_r.get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)
    j_coupling = r_r.lumped_model_get_j(C1, L1, C2+Cs, L2, Cg, Lg)
    #omegas = np.linspace(10*np.pi*1e9, 24*np.pi*1e9, 10000)
    #Zs = r_r.lumped_model_Z21(omegas, C1, L1, C2+Cs, L2, Cg, Lg)
    #Zs = np.abs(Zs)


    total_error = (f1-target[0])**2 + (1e13*(len_r2-target[1]))**2 + (5*(notch_f-target[2]))**2 + (5e1*(j_coupling-target[3]))**2
    #print("f1 = ", f1/1e9)
    #print("len_r2 = ", len_r2*1e6)
    #print("target[2] = ", target[1]*1e6 )
    return total_error

def print_progress1(xk, convergence):
    #print("xk = ", xk)
    print("\r convergence = ", "{:.5f}".format(convergence), "  ", " Lgn =",int(xk[0]*1e7)/10, " Lgf =",int(xk[1]*1e7)/10, " Lc =",int(xk[2]*1e7)/10, " Lrn =",int(xk[3]*1e7)/10, " Lrf =",int(xk[4]*1e7)/10, " d =",int(xk[5]*1e7)/10, "     ",end="")
    sys.stdout.flush()

def print_progress2(xk):
    convergence=0
    #print("xk = ", xk)
    print("\r convergence = ", "{:.5f}".format(convergence), "  ", " Lgn =",int(xk[0]*1e7)/10, " Lgf =",int(xk[1]*1e7)/10, " Lc =",int(xk[2]*1e7)/10, " Lrn =",int(xk[3]*1e7)/10, " Lrf =",int(xk[4]*1e7)/10, " d =",int(xk[5]*1e7)/10, "     ",end="")
    sys.stdout.flush()

def print_message(sol, target, Cs, calibration_len, plot=0):
    solution = sol.x
    result = function(solution, [0,0,0,0], Cs, calibration_len, weighted_error=0)
    
    print("\n \n")
    print("################################")
    print("    r_r: Solution found")
    print(" Lgn = ", int(solution[0]*1e7)/10, "[um]")
    print(" Lgf = ", int(solution[1]*1e7)/10, "[um]")
    print(" Lc = ", int(solution[2]*1e7)/10, "[um]")
    print(" Lrn = ", int(solution[3]*1e7)/10, "[um]")
    print(" Lrf = ", int(solution[4]*1e7)/10, "[um]")
    print(" d = ", int(solution[5]*1e7)/10, "[um]")
    print("-------------------------------")
    print("         Results in")
    print(" fr1 = ", int(result[0]*1e-6)/1000, "[GHz]")
    print(" len_r2 = ", int(result[1]*1e7)/10, "[um]")
    print(" fn = ", int(result[2]*1e-6)/1000, "[GHz]")
    print(" J = ", int(result[3]*1e-3/2/np.pi)/1000, "[MHz]")
    print("################################")

    if plot==1:
        omegas = np.linspace(10*np.pi*1e9, 24*np.pi*1e9, 10000)
        C1, L1, C2, L2, Cg, Lg = r_r.get_lumped_elements(solution[2], solution[1]+calibration_len/2, solution[0]+calibration_len/2, solution[4]+calibration_len/2, solution[3]+calibration_len/2, cap.get_Lm(solution[5]), cap.get_Cm(solution[5]))
        Zs = r_r.lumped_model_Z21(omegas, C1, L1, C2+Cs, L2, Cg, Lg)
        Zs = np.abs(Zs)
        
        plt.plot(omegas/(2*np.pi*1e9), Zs, label="lumped element model", color="blue")
        plt.vlines([target[0]*1e-9,cf.lambda_by_4_f(target[1])*1e-9], 1e-3, 1e5, color="orange", label="resonators")
        plt.vlines([target[2]*1e-9], 1e-3, 1e5, color="red", label="notch")
        plt.yscale("log")
        plt.show()
    return 0

def solve_for_r_r(target, x0, Cs, calibration_len, show=0, plot=0):
    print("\n")
    print("Started solving for r_r...")
    #sol = opt.root(function, x0, args=(target, Cs))
    bounds=[(750e-6,1500e-6),(10e-6,2500e-6),(200e-6,400e-6),(750e-6,1500e-6),(10e-6,2500e-6),(3e-6,10e-6)]
    #
    sol = opt.differential_evolution(function2, args=(target, Cs, calibration_len), bounds=bounds, callback=print_progress1)
    #sol = opt.minimize(function2, args=(target, Cs, calibration_len), x0=x0, bounds=bounds, callback=print_progress2)
    if show!=0:
        print_message(sol, target, Cs, calibration_len, plot=plot)
    return sol


def get_notch_frequency_derivative(l_c, l_Gf, l_Gn, l_Rf, l_Rn, d, jump_length=50e-6):
    # NOTE: calibration needs to be implemented in this too !!!
    Cm_per_len = cap.get_Cm(d)
    Lm_per_len = cap.get_Lm(d)
    upper_notch_f = r_r.find_notch_filter_frequency(l_c+jump_length, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
    lower_notch_f = r_r.find_notch_filter_frequency(l_c-jump_length, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
    print("Derivative by lc: ", (upper_notch_f-lower_notch_f)/1e6)
    
    upper_notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf+jump_length, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
    lower_notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf-jump_length, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
    print("Derivative by lgf: ", (upper_notch_f-lower_notch_f)/1e6)
    
    upper_notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn+jump_length, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
    lower_notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn-jump_length, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
    print("Derivative by lgn: ", (upper_notch_f-lower_notch_f)/1e6)
    
    upper_notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf+jump_length, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
    lower_notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf-jump_length, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
    print("Derivative by lrf: ", (upper_notch_f-lower_notch_f)/1e6)
    
    upper_notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn+jump_length, Lm_per_len, Cm_per_len)/(2*np.pi)
    lower_notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn-jump_length, Lm_per_len, Cm_per_len)/(2*np.pi)
    print("Derivative by lrn: ", (upper_notch_f-lower_notch_f)/1e6)
    
    Cm_per_len = cap.get_Cm(d+0.5e-6)
    Lm_per_len = cap.get_Lm(d+0.5e-6)
    upper_notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
    Cm_per_len = cap.get_Cm(d-0.5e-6)
    Lm_per_len = cap.get_Lm(d-0.5e-6)
    lower_notch_f = r_r.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)/(2*np.pi)
    print("Derivative by d: ", (upper_notch_f-lower_notch_f)/1e6)
    
    return 0

def plot_solution_space():
    num=100
    bounds=[(780e-6,2500e-6),(10e-6,2500e-6),(10e-6,3000e-6),(770e-6,2500e-6),(10e-6,2500e-6),(2e-6,50e-6)]
    
    range_l_Gn = np.linspace(bounds[0][0],bounds[0][1], num)
    range_l_Gf = np.linspace(bounds[1][0],bounds[1][1], num)
    range_l_c = np.linspace(bounds[2][0],bounds[2][1], num)
    range_l_Rn = np.linspace(bounds[3][0],bounds[3][1], num)
    range_l_Rf = np.linspace(bounds[4][0],bounds[4][1], num)
    range_d = np.linspace(bounds[5][0],bounds[5][1], num)

    Lgn =  790.0e-6
    Lgf =  1566.1e-6
    Lc =  416.3e-6
    Lrn =  790.0e-6
    Lrf =  1341.6e-6
    d =  3.0e-6
    Cm_per_len = cap.get_Cm(d)
    Lm_per_len = cap.get_Lm(d)
     
    values = np.zeros((num,num))
    inputs = [0,0,0,0,0,0]
    target = [0,0,8.5e9,0]
    calibration_len=250e-6
    Cs = 0
    for i in range(num):
        for j in range(num):
            #print(inputs)
            inputs[0] = Lgn+calibration_len/2
            inputs[1] = range_l_Gf[i]+calibration_len/2
            inputs[2] = Lc
            inputs[3] = Lrn+calibration_len/2
            inputs[4] = range_l_Rf[j]+calibration_len/2
            inputs[5] = d
            values[i,j] = function(inputs, target, Cs, calibration_len)[2]
    
    plt.pcolormesh(range_l_Gf*1e6, range_l_Rf*1e6, values)
    plt.colorbar()
    plt.show()
    print(values)


#plot_solution_space()
"""
x0 = [800e-6, 1600e-6, 370e-6, 827e-6, 1100e-6, 5e-6]
target=[10e9, 2863e-6, 8e9, 50e6*2*np.pi] # target format: Fr1, len_res2, Fn, J
solve_for_r_r(target,x0,1,1)
"""



