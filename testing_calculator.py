
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.optimize as opt
import scipy.interpolate as interpol

import common_formulas as cf
import q_r_formulas as qr_f
import q_r_algo as qr_a
import r_r_formulas as rr_f
import r_r_algo as rr_a
import cap_util as cap




def plot_transmission(l_Rf, l_Rn, l_Gf, l_Gn, l_c, d):

    omegas = np.linspace(5*2*np.pi*1e9, 12*2*np.pi*1e9, 10000)
    
    Cm_per_len = cap.get_Cm(d)
    Lm_per_len = cap.get_Lm(d)

    C1, L1, C2, L2, Cg, Lg = rr_f.get_lumped_elements(Lc, Lgf, Lgn, Lrf, Lrn, Lm_per_len, Cm_per_len)
    Zs = rr_f.lumped_model_transmission(C1, L1, C2, L2, Cg, Lg, omegas)
    Zs = np.abs(Zs)
    
    # notches
    notch_freq = rr_f.find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm_per_len, Cm_per_len)
    print("notch_freq = " + str(int(notch_freq/1e6/2/np.pi)/1e3))
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


    
#
# Qubit frequency
#

def get_Ec(C):
    e = 1.6022e-19
    return e**2/2/C

def get_Ej(L):
    phi0 = 2.0678e-15
    return ((phi0/2/np.pi)**2)/L

def omega_q(L, C):
    hbar = 1.05457182e-34
    Ec = get_Ec(C)
    Ej = get_Ej(L)
    epsilon = np.sqrt((2*Ec/Ej))
    om = np.sqrt(8*Ec*Ej) - Ec*(1 + 0.25*epsilon + (21/(2**7))*(epsilon**2) + (19/(2**7))*(epsilon**3) + (5319/(2**15))*epsilon**4)
    return om/hbar

#print(omega_q(5.5e-9, 65e-15)/2/np.pi/1e9)

def chi(g, Ec, fq, fr):
    hbar = 1.05457182e-34
    om_q = fq *2 *np.pi
    om_r = fr *2 *np.pi
    d = -om_r + om_q
    return -g**2 * Ec/hbar/(d*(d-Ec/hbar))

def chi_opt(g, Ec, fq, fr, target):
    return chi(g, Ec, fq, fr) - target
    
Cq = 65e-15
g0 = 350e6*2*np.pi
fq = 8.1e9
fr = 10.2e9

Cqs = np.array([65, 53.5, 55, 66.5])*1e-15
fqs = np.array([8.1, 8.9, 8.8, 8])*1e9
frs = np.array([10.16, 10.5, 10.33, 10])*1e9


chi_target = -10*np.pi*1e6
"""
for i in range(4):
    print("---------------")
    print(i)
    Ec = get_Ec(Cqs[i])
    g = opt.fsolve(chi_opt, x0=g0, args=(Ec, fqs[i], frs[i], chi_target))
    print("g/2/pi = ", g/2/np.pi/1e6)
    print("2chi/2/pi = ", chi(g,Ec,fqs[i],frs[i])/np.pi/1e6)

    r = qr_f.g2r(g, fqs[i]*2*np.pi, frs[i]*2*np.pi, Cqs[i])
    print("r/2/pi = ", r/1e5)

quit()

"""
"""
#
# J coupling
#

# manual formula for converting the r to g
g = np.array([350, 255, 250, 345])*1e6*2*np.pi
om_q = np.array([8.1, 8.9, 8.8, 8])*1e9*2*np.pi
om_r = np.array([10.2, 10.5, 10.4, 10.1])*1e9*2*np.pi
c_q = np.array([65, 53.5, 54.75, 66.5])*1e-15


for i in range(4):
    print(qr_f.g2r(g[i], om_q[i], om_r[i], c_q[i])/1e5)
quit()
"""


"""

# manual formula for converting the r to g

# f1 r1 f2 r2
#simulations = np.array([[9.68e9, 35028, 9.77e9, 1.472e5],[9.54e9, 86448, 9.578e9, 1.202e5],[9.508e9, 12885, 9.7319e9, 1.496e5],[9.46e9, 22271, 9.58e9, 1.475e5]])
#simulations = np.array([[9.62e9, 15534, 9.9142e9, 1.53e5],[9.61e9, 9829, 9.965e9, 1.553e5],[9.608e9, 9829.3, 9.965e9, 1.553e5],[9.617e9,9787,9.9925e9,1.53e5],[9.604e9,10991, 9.915e9, 1.547e5], [9.6232e9, 12302, 9.442e9, 1.6314e5],[9.6218e9,13019,9.989e9,1.533e5],[9.568e9, 11578,9.97e9,1.577e5],[9.569e9,10885,9.9537e9,1.5894e5]])

#simulations = np.array([[9.512e9, 8299.7, 9.8004e9, 1.4599e5],[9.5252e9, 8555.7, 9.7976e9, 1.4726e5], [9.5642e9, 7085.6, 9.9076e9, 1.503e5], [9.629e9, 12424, 9.8995e9, 1.5129e5]])

#simulations = np.array([[9.938e9,54817,9.956e9,1.4e5],[9.8655e9,43422,9.9221e9,1.467e5],[9.787e9,57736,9.7943e9,1.415e5],[9.7174e9,53302,9.689e9,1.3982e5]])
#simulations = np.array([[9.925e9,4873,9.232e9,1.4667e5],[9.854e9,5040,9.1913e9,1.47e5],[9.8e9,4621,9.122e9,1.43e5],[9.711e9,4993,9.0334e9,1.43e5]])
#simulations = np.array([[9.836e9, 2217, 8.85e9, 1.408e5],[9.844e9, 4016, 8.864e9, 1.464e5],[9.9e9, 5110, 8.896e9, 1.46e5]])
#simulations = np.array([[9.836e9, 2217, 8.9343e9, 1.4323e5],[9.844e9, 4016, 8.864e9, 1.464e5],[9.9e9, 5110, 8.896e9, 1.46e5],[9.91e9, 5893, 8.9e9,1.64e5]])

"""
"""
J_couplings = np.array([0.02200551, 0.02957635, 0.03718251, 0.03710093])*1e3
names = ["fine", "finer", "extra fine", "extremely fine"]
plt.plot(names, J_couplings, linestyle="none", marker="o")
plt.ylabel("J_coupling [MHz]")
plt.show()
quit()
"""


def linear(x, a, b):
    return x*a+b

def fitfunc(x, middle, minim, a1, a2):
    if x<middle:
        return (middle-x)*a1+minim
    if x>middle:
        return (x-middle)*a2+minim
    if x==middle:
        return minim

def ev_func(ind, J, om1, a, b):
    om2 = a*ind+b
    lambda1 = 0.5*(om1 + om2 - np.sqrt((om1+om2)**2 - 4*(om1*om2-J**2)))
    lambda2 = 0.5*(om1 + om2 + np.sqrt((om1+om2)**2 - 4*(om1*om2-J**2)))
    return np.array([lambda1, lambda2])

def residual(params, ind, lambdas):
    J = params[0]
    om1 = params[1]
    a = params[2]
    b = params[3]

    errors = np.zeros(len(ind))
    
    for i in range(len(ind)):
        [l1,l2] = ev_func(ind[i], J, om1, a, b)
        errors[i] = (lambdas[i,0]-l1)**2 + (lambdas[i,1]-l2)**2
    return errors
        
#filename = "C:\\Users\\mluka\\work\\data\\j_coupling_mesh_fine.csv"
filename = "C:\\Users\\mluka\\work\\data\\transfer.csv"
data = np.genfromtxt(filename,delimiter="," ,skip_header=5)
#data = data[16:-20]

#filename = "C:\\Users\\mluka\\work\\data\\j_coupling_mesh_finer.csv"
#filename = "C:\\Users\\mluka\\work\\data\\transfer.csv"
#data2 = np.genfromtxt(filename,delimiter="," ,skip_header=5)

#filename = "C:\\Users\\mluka\\work\\data\\j_coupling_mesh_extrafine.csv"
#filename = "C:\\Users\\mluka\\work\\data\\transfer.csv"
#data3 = np.genfromtxt(filename,delimiter="," ,skip_header=5)

#filename = "C:\\Users\\mluka\\work\\data\\j_coupling_mesh_extraextra.csv"
#filename = "C:\\Users\\mluka\\work\\data\\transfer.csv"
#data4 = np.genfromtxt(filename,delimiter="," ,skip_header=5)

def get_J_method1(data):
    inds = data[::2,0]
    lambdas = np.zeros((len(inds),2))
    lambdas[:,0] = data[::2,2]
    lambdas[:,1] = data[1::2,2]
    f1 = interpol.interp1d(inds, lambdas[:,0], "quadratic", fill_value="extrapolate")
    f2 = interpol.interp1d(inds, lambdas[:,1], "quadratic", fill_value="extrapolate")

    def f3(ind):
        return f2(ind)-f1(ind)
    xmin = opt.minimize(f3, x0=0.04)
    J = f3(xmin.x)/2
    return J, xmin.x, (f1(xmin.x)+f2(xmin.x))/2

def get_J_method2(data):
    inds = data[::2,0]
    lambdas = np.zeros((len(inds),2))
    lambdas[:,0] = data[::2,2]
    lambdas[:,1] = data[1::2,2]

    result = opt.least_squares(residual,x0=[0.03,10,-10,10],args=(inds, lambdas))
    return(result.x)

J_fine, x_fine, y_fine = get_J_method1(data)
#J_finer, x_finer, y_finer = get_J_method1(data2)
#J_extrafine, x_extrafine, y_extrafine = get_J_method1(data3)
#J_extraextra, x_extraextra, y_extraextra = get_J_method1(data4)
print(J_fine)
#print(J_finer)
#print(J_extrafine)
#print(J_extraextra)
J_fine = get_J_method2(data)
#J_finer = get_J_method2(data2)
#J_extrafine = get_J_method2(data3)
#J_extraextra = get_J_method2(data4)
print(J_fine)
#print(J_finer)
#print(J_extrafine)
#print(J_extraextra)

#inds_continous = np.linspace(inds[0],inds[-1],1000)
#lambs = ev_func(inds_continous, result.x[0], result.x[1], result.x[2], result.x[3])
#lambs = ev_func(inds_continous, 0, om1, a, b)
#print(lambs)


#plt.plot(data[::2,0],differences,linestyle="none",marker="o")
#plt.plot(ds_left,fitfunc0(ds_left, *left_popt))
#plt.plot(ds_right,fitfunc0(ds_right, *right_popt))
#plt.plot(x_min,y_min,linestyle="none",marker="*",color="red")
#plt.plot(ds, np.vectorize(fitfunc)(ds,*popt))
plt.plot(data[:,0],data[:,2], linestyle="none", marker="o", label="fine")
#plt.plot(data2[:,0],data2[:,2], linestyle="none", marker="o", label="finer")
#plt.plot(data3[:,0],data3[:,2], linestyle="none", marker="o", label="extrafine")
#plt.plot(data4[:,0],data4[:,2], linestyle="none", marker="o", label="extraextra")
#plt.plot(inds_continous, linear(inds_continous, result.x[2], result.x[3]))
#plt.plot(inds_continous,lambs[0])
#plt.plot(inds_continous,lambs[1])
#plt.plot(inds_continous,f1(inds_continous))
#plt.plot(inds_continous,f2(inds_continous))
plt.plot(x_fine, y_fine, linestyle="none", marker="x", color="red")
#plt.plot(x_finer, y_finer, linestyle="none", marker="x", color="red")
#plt.plot(x_extrafine, y_extrafine, linestyle="none", marker="x", color="red")
#plt.plot(x_extraextra, y_extraextra, linestyle="none", marker="x", color="red")
plt.xlabel("Lumped inductance [nH]")
plt.ylabel("Eigenfrequencies [GHz]")
plt.legend()
plt.show()

#plt.plot([20, 25, 30],[23, 26.4, 29.5], linestyle="none", marker="o", label="interpolation")
#plt.plot([20, 25, 30],[23, 24, 29.9], linestyle="none", marker="v", label="fit")
#plt.plot(np.linspace(20,30,100), np.linspace(20,30,100), linestyle="--", color="grey")
#plt.xlabel("Targeted J [MHz]")
#plt.ylabel("Simulated J [MHz]")
#plt.show()

"""
om_q = 8.1e9*2*np.pi
c_q = 7.5e-14

J_arr = np.zeros(int(len(data)/2))

for i in range(int(len(data)/2)):
    print("----------")

    r1 = data[2*i,2]
    om_1 = data[2*i,1]*2*np.pi*1e9
    g1 = qr_f.r2g(r1, om_q, om_1, c_q)
    print(g1/2/np.pi/1e6)
    
    r2 = data[2*i+1,2]
    om_2 = data[2*i+1,1]*2*np.pi*1e9
    g2 = qr_f.r2g(r2, om_q, om_2, c_q)
    print(g2/2/np.pi/1e6)

    delta1 = om_q-om_1
    delta2 = om_2-om_1
    
    J = (g2/g1)*(2*delta1*delta2/(delta1+delta2))
    J_arr[i] = -J/2/np.pi/1e6
    
    print(J/2/np.pi/1e6)




"""

# Lumped element model
calibration_len=400e-6
"""
#peters q1
Lgn =  1192e-6 +calibration_len/2
Lgf =  1338e-6 +calibration_len/2
Lc =  338e-6
Lrn =  908e-6 +calibration_len/2
Lrf =  1508e-6 +calibration_len/2
d = 5e-6
"""

# # peters q2
# Lgn =  1212e-6 +calibration_len/2
# Lgf =  1279e-6 +calibration_len/2
# Lc =  338e-6
# Lrn =  926e-6 +calibration_len/2
# Lrf =  1418e-6 +calibration_len/2
# d = 5e-6



# peters q4
# Lgn =  1199e-6 +calibration_len/2
# Lgf =  1383e-6 +calibration_len/2
# Lc =  338e-6
# Lrn =  908e-6 +calibration_len/2
# Lrf =  1593e-6 +calibration_len/2
# d = 5e-6



# # peters q3
#Lgn =  1203e-6 +calibration_len/2
#Lgf =  1630e-6 +calibration_len/2
#Lc =  338e-6
#Lrn =  920e-6 +calibration_len/2
#Lrf =  1468e-6 +calibration_len/2
#d = 5e-6


# #q1
# Lgn =  770e-6 +calibration_len/2
# Lgf =  1420e-6 +calibration_len/2
# Lc =  304e-6
# Lrn =  665e-6 +calibration_len/2
# Lrf =  1390e-6 +calibration_len/2
# d = 7e-6

# #q1 mux2
# Lgn =  720e-6 +calibration_len/2
# Lgf =  1460e-6 +calibration_len/2
# Lc =  304e-6
# Lrn =  625e-6 +calibration_len/2
# Lrf =  1430e-6 +calibration_len/2
# d = 5e-6

# q2
Lgn =  780e-6 +calibration_len/2
Lgf =  1325e-6 +calibration_len/2
Lc =  363e-6
Lrn =  665e-6 +calibration_len/2
Lrf =  1290e-6 +calibration_len/2
d = 7e-6

# # mux 2 
# # q2
# Lgn =  1015e-6 +calibration_len/2
# Lgf =  1090e-6 +calibration_len/2
# Lc =  363e-6
# Lrn =  880e-6 +calibration_len/2
# Lrf =  1076e-6 +calibration_len/2
# d = 3e-6




# # q3
# Lgn =  850e-6 +calibration_len/2
# Lgf =  1330e-6 +calibration_len/2
# Lc =  344e-6
# Lrn =  710e-6 +calibration_len/2
# Lrf =  1322e-6 +calibration_len/2
# d = 5e-6

# #mux 2
# # q3
# Lgn =  950e-6 +calibration_len/2
# Lgf =  1230e-6 +calibration_len/2
# Lc =  344e-6
# Lrn =  810e-6 +calibration_len/2
# Lrf =  1222e-6 +calibration_len/2
# d = 3e-6

# # q4
# Lgn =  790e-6 +calibration_len/2
# Lgf =  1470e-6 +calibration_len/2
# Lc =  294e-6
# Lrn =  685e-6 +calibration_len/2
# Lrf =  1360e-6 +calibration_len/2
# d = 8e-6

# #mux 2
# # q4
# Lgn =  690e-6 +calibration_len/2
# Lgf =  1570e-6 +calibration_len/2
# Lc =  294e-6
# Lrn =  625e-6 +calibration_len/2
# Lrf =  1420e-6 +calibration_len/2
# d = 5e-6


Cs = 2.72e-14
#Cs = 0
#Rs = 10206.01
#Lextra = 0.2e-9
#Lextra=0



Cm_per_len = cap.get_Cm(d)
Lm_per_len = cap.get_Lm(d)

for shift in [0]:

    C1, L1, C2, L2, Cg, Lg = rr_f.get_lumped_elements(Lc, Lgf-shift, Lgn+shift, Lrf+shift, Lrn-shift, Lm_per_len, Cm_per_len)
    # print("C1 = ", C1)
    # print("L1 = ", L1)
    # print("C2 = ", C2)
    # print("L2 = ", L2)
    # print("Cg = ", Cg)
    # print("Lg = ", Lg)
    print("shift: ", shift*1e6)
    notch_freq = rr_f.find_notch_filter_frequency(Lc, Lgf-shift, Lgn+shift, Lrf+shift, Lrn-shift, Lm_per_len, Cm_per_len)
    print("notch_freq = " + str(int(notch_freq/1e6/2/np.pi)/1e3))
    j_coupling = rr_f.lumped_model_get_j(C1, L1, C2+Cs, L2, Cg, Lg)
    #plot_transmission(Lrf, Lrn, Lgf, Lgn, Lc, d)
    print(j_coupling/2/np.pi/1e6)




#print(j_coupling/1e6/np.pi/2)

#Zs = lumped_model_Z21(omegas, C1, L1, C2+Cs, L2, Cg, Lg)
#Zs = np.abs(Zs)
        
#plt.plot(omegas/(2*np.pi*1e9), Zs, label="lumped element model", color="blue")
#plt.yscale("log")
#plt.show()

    
#plt.plot(data[::2,0],J_arr/8.82, linestyle="none", marker="o")
#plt.plot(Lc*1e6,J_lumped)
#plt.xlabel("Lc [um]")
#plt.ylabel("J/2pi [MHz]")
#plt.show()


#
# Notch frequency derivative stuff
#


"""
jump_lengths = np.linspace(0,3e-6,100)
upper_freqs = np.zeros(len(jump_lengths))
lower_freqs = np.zeros(len(jump_lengths))

calibration = 470e-6

Lgn =  860.0e-6
Lgf =  1278.8e-6
Lc =  457.7e-6
Lrn =  860.5e-6
Lrf =  1586.8e-6
d=5e-6

get_notch_fequency_derivative(Lc, Lgf+calibration/2, Lgn+calibration/2, Lrf+calibration/2, Lrn+calibration/2, d, 50e-6)
"""

