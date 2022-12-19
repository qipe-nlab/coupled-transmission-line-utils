# This file takes care of fitting the values obtained from the 2d simulations

######################
### IMPORTANT NOTE ###
######################

# The functions that should be used in other scripts are get_Cm and get_Lm

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as lin
import scipy.interpolate as interpol


def generate_Cm(sample_name="capacitance_d-1-50-100_l-1000.csv"):

    filename = "C:\\Users\\mluka\\work\\data\\capind\\" + sample_name
    data = np.genfromtxt(filename,delimiter="," ,skip_header=5)
    len_data = int(np.shape(data)[0]/2)

    # Format: [[d,Cm],[d,Cm],[d,Cm]]
    Cm_mat = np.zeros((len_data, 2))

    for i in range(len_data):
        Cm_mat[i,0] = data[2*i,0]
        Cm_mat[i,1] = data[2*i,2]/1e-3 # account for simulation being of 1mm

    return Cm_mat


def generate_Lm(sample_name="mcapacitance_d-1-50-100_l-1000.csv"):

    epsilon = 8.854e-12
    mu = 1.256e-6
    
    filename = "C:\\Users\\mluka\\work\\data\\capind\\" + sample_name
    data = np.genfromtxt(filename,delimiter="," ,skip_header=5)
    len_data = int(np.shape(data)[0]/2)

    # Format: [[d,Lm],[d,Lm],[d,Lm]]
    Lm_mat = np.zeros((len_data, 2))

    for i in range(len_data):
        C_matrix = np.array([[data[2*i,1], data[2*i,2]],[data[2*i+1,1], data[2*i+1,2]]])/1e-3
        L_matrix = lin.inv(C_matrix)*epsilon*mu

        Lm_mat[i,0] = data[2*i,0]
        Lm_mat[i,1] = L_matrix[0,1]

    return Lm_mat


def get_Cm_nonvec(d):

    Cm_mat = np.array([[1.00000000e+00, 1.50370318e-11], [2.00000000e+00, 1.28176495e-11], [3.00000000e+00, 1.12177246e-11], [4.00000000e+00, 9.96369927e-12], [5.00000000e+00, 8.94052811e-12], [6.00000000e+00, 8.08682762e-12], [7.00000000e+00, 7.36199016e-12], [8.00000000e+00, 6.73831905e-12], [9.00000000e+00, 6.19573028e-12], [1.00000000e+01, 5.72086080e-12], [1.10000000e+01, 5.30132999e-12], [1.20000000e+01, 4.92745225e-12], [1.30000000e+01, 4.59423231e-12], [1.40000000e+01, 4.29470251e-12], [1.50000000e+01, 4.02450889e-12], [1.60000000e+01, 3.77978074e-12], [1.70000000e+01, 3.55708594e-12], [1.80000000e+01, 3.35442782e-12], [1.90000000e+01, 3.16833752e-12], [2.00000000e+01, 2.99802731e-12]])

    if d<Cm_mat[0,0]:
        raise ValueError('No prediction for Cm, d is too small.')
    if d>=Cm_mat[-1,0]:
        raise ValueError('No prediction for Cm, d is too small.')

    # find between which two datapoints our requested d is
    index = 0
    while d>Cm_mat[index+1, 0]:
        index += 1

    # linear interpolation between points
    dx = Cm_mat[index+1, 0] - Cm_mat[index, 0]
    dydx = (Cm_mat[index+1, 1] - Cm_mat[index, 1])/dx
    
    Cm = (Cm_mat[index, 1]) + dydx*(d-Cm_mat[index,0])

    return Cm

def get_Cm_interpol(d):
    Cm_mat = np.array([[1.00000000e+00, 1.50370318e-11], [1.50000000e+00, 1.38187520e-11], [2.00000000e+00, 1.28176495e-11], [2.50000000e+00, 1.19642927e-11], [3.00000000e+00, 1.12177246e-11], [3.50000000e+00, 1.05558056e-11], [4.00000000e+00, 9.96369927e-12], [4.50000000e+00, 9.42783813e-12], [5.00000000e+00, 8.94052811e-12], [5.50000000e+00, 8.49606695e-12], [6.00000000e+00, 8.08682762e-12], [6.50000000e+00, 7.71056795e-12], [7.00000000e+00, 7.36199016e-12], [7.50000000e+00, 7.03856503e-12], [8.00000000e+00, 6.73831905e-12], [8.50000000e+00, 6.45844624e-12], [9.00000000e+00, 6.19573028e-12], [9.50000000e+00, 5.95092974e-12], [1.00000000e+01, 5.72086080e-12], [1.05000000e+01, 5.50431223e-12], [1.10000000e+01, 5.30132999e-12], [1.15000000e+01, 5.10927178e-12], [1.20000000e+01, 4.92745225e-12], [1.25000000e+01, 4.75662738e-12], [1.30000000e+01, 4.59423231e-12], [1.35000000e+01, 4.44082587e-12], [1.40000000e+01, 4.29470251e-12], [1.45000000e+01, 4.15623864e-12], [1.50000000e+01, 4.02450889e-12], [1.55000000e+01, 3.89957092e-12], [1.60000000e+01, 3.77978074e-12], [1.65000000e+01, 3.66599253e-12], [1.70000000e+01, 3.55708594e-12], [1.75000000e+01, 3.45341660e-12], [1.80000000e+01, 3.35442782e-12], [1.85000000e+01, 3.25928889e-12], [1.90000000e+01, 3.16833752e-12], [1.95000000e+01, 3.08118740e-12], [2.00000000e+01, 2.99802731e-12], [2.05000000e+01, 2.91774255e-12], [2.10000000e+01, 2.84102594e-12], [2.15000000e+01, 2.76720261e-12], [2.20000000e+01, 2.69623703e-12], [2.25000000e+01, 2.62806109e-12], [2.30000000e+01, 2.56226059e-12], [2.35000000e+01, 2.49936522e-12], [2.40000000e+01, 2.43825736e-12], [2.45000000e+01, 2.37982800e-12], [2.50000000e+01, 2.32319135e-12], [2.55000000e+01, 2.26879000e-12], [2.60000000e+01, 2.21617328e-12], [2.65000000e+01, 2.16539854e-12], [2.70000000e+01, 2.11616745e-12], [2.75000000e+01, 2.06897809e-12], [2.80000000e+01, 2.02305598e-12], [2.85000000e+01, 1.97873162e-12], [2.90000000e+01, 1.93570311e-12], [2.95000000e+01, 1.89433481e-12], [3.00000000e+01, 1.85407274e-12], [3.05000000e+01, 1.81510429e-12], [3.10000000e+01, 1.77739355e-12], [3.15000000e+01, 1.74094157e-12], [3.20000000e+01, 1.70534396e-12], [3.25000000e+01, 1.67098021e-12], [3.30000000e+01, 1.63767627e-12], [3.35000000e+01, 1.60527678e-12], [3.40000000e+01, 1.57375707e-12], [3.45000000e+01, 1.54342609e-12], [3.50000000e+01, 1.51378668e-12], [3.55000000e+01, 1.48485275e-12], [3.60000000e+01, 1.45685660e-12], [3.65000000e+01, 1.42960162e-12], [3.70000000e+01, 1.40309872e-12], [3.75000000e+01, 1.37724135e-12], [3.80000000e+01, 1.35229549e-12], [3.85000000e+01, 1.32786810e-12], [3.90000000e+01, 1.30421364e-12], [3.95000000e+01, 1.28089509e-12], [4.00000000e+01, 1.25853942e-12], [4.05000000e+01, 1.23652075e-12], [4.10000000e+01, 1.21502824e-12], [4.15000000e+01, 1.19449342e-12], [4.20000000e+01, 1.17400121e-12], [4.25000000e+01, 1.15413225e-12], [4.30000000e+01, 1.13478322e-12], [4.35000000e+01, 1.11598960e-12], [4.40000000e+01, 1.09764123e-12], [4.45000000e+01, 1.07962483e-12], [4.50000000e+01, 1.06202586e-12], [4.55000000e+01, 1.04493374e-12], [4.60000000e+01, 1.02821495e-12], [4.65000000e+01, 1.01187883e-12], [4.70000000e+01, 9.96036092e-13], [4.75000000e+01, 9.80360455e-13], [4.80000000e+01, 9.65162058e-13], [4.85000000e+01, 9.50366563e-13], [4.90000000e+01, 9.35695396e-13], [4.95000000e+01, 9.21522663e-13], [5.00000000e+01, 9.07535237e-13]])
    f = interpol.interp1d(Cm_mat[:,0]/1e6, Cm_mat[:,1], "quadratic")
    return f(d)

def get_Cm(d):
    return np.vectorize(get_Cm_interpol)(d)


def get_Cm_interpol_lin(d):
    # Dont use this function, this was just used for comparing linear and quadratic interpolation
    Cm_mat = np.array([[1.00000000e+00, 1.50370318e-11], [1.50000000e+00, 1.38187520e-11], [2.00000000e+00, 1.28176495e-11], [2.50000000e+00, 1.19642927e-11], [3.00000000e+00, 1.12177246e-11], [3.50000000e+00, 1.05558056e-11], [4.00000000e+00, 9.96369927e-12], [4.50000000e+00, 9.42783813e-12], [5.00000000e+00, 8.94052811e-12], [5.50000000e+00, 8.49606695e-12], [6.00000000e+00, 8.08682762e-12], [6.50000000e+00, 7.71056795e-12], [7.00000000e+00, 7.36199016e-12], [7.50000000e+00, 7.03856503e-12], [8.00000000e+00, 6.73831905e-12], [8.50000000e+00, 6.45844624e-12], [9.00000000e+00, 6.19573028e-12], [9.50000000e+00, 5.95092974e-12], [1.00000000e+01, 5.72086080e-12], [1.05000000e+01, 5.50431223e-12], [1.10000000e+01, 5.30132999e-12], [1.15000000e+01, 5.10927178e-12], [1.20000000e+01, 4.92745225e-12], [1.25000000e+01, 4.75662738e-12], [1.30000000e+01, 4.59423231e-12], [1.35000000e+01, 4.44082587e-12], [1.40000000e+01, 4.29470251e-12], [1.45000000e+01, 4.15623864e-12], [1.50000000e+01, 4.02450889e-12], [1.55000000e+01, 3.89957092e-12], [1.60000000e+01, 3.77978074e-12], [1.65000000e+01, 3.66599253e-12], [1.70000000e+01, 3.55708594e-12], [1.75000000e+01, 3.45341660e-12], [1.80000000e+01, 3.35442782e-12], [1.85000000e+01, 3.25928889e-12], [1.90000000e+01, 3.16833752e-12], [1.95000000e+01, 3.08118740e-12], [2.00000000e+01, 2.99802731e-12], [2.05000000e+01, 2.91774255e-12], [2.10000000e+01, 2.84102594e-12], [2.15000000e+01, 2.76720261e-12], [2.20000000e+01, 2.69623703e-12], [2.25000000e+01, 2.62806109e-12], [2.30000000e+01, 2.56226059e-12], [2.35000000e+01, 2.49936522e-12], [2.40000000e+01, 2.43825736e-12], [2.45000000e+01, 2.37982800e-12], [2.50000000e+01, 2.32319135e-12], [2.55000000e+01, 2.26879000e-12], [2.60000000e+01, 2.21617328e-12], [2.65000000e+01, 2.16539854e-12], [2.70000000e+01, 2.11616745e-12], [2.75000000e+01, 2.06897809e-12], [2.80000000e+01, 2.02305598e-12], [2.85000000e+01, 1.97873162e-12], [2.90000000e+01, 1.93570311e-12], [2.95000000e+01, 1.89433481e-12], [3.00000000e+01, 1.85407274e-12], [3.05000000e+01, 1.81510429e-12], [3.10000000e+01, 1.77739355e-12], [3.15000000e+01, 1.74094157e-12], [3.20000000e+01, 1.70534396e-12], [3.25000000e+01, 1.67098021e-12], [3.30000000e+01, 1.63767627e-12], [3.35000000e+01, 1.60527678e-12], [3.40000000e+01, 1.57375707e-12], [3.45000000e+01, 1.54342609e-12], [3.50000000e+01, 1.51378668e-12], [3.55000000e+01, 1.48485275e-12], [3.60000000e+01, 1.45685660e-12], [3.65000000e+01, 1.42960162e-12], [3.70000000e+01, 1.40309872e-12], [3.75000000e+01, 1.37724135e-12], [3.80000000e+01, 1.35229549e-12], [3.85000000e+01, 1.32786810e-12], [3.90000000e+01, 1.30421364e-12], [3.95000000e+01, 1.28089509e-12], [4.00000000e+01, 1.25853942e-12], [4.05000000e+01, 1.23652075e-12], [4.10000000e+01, 1.21502824e-12], [4.15000000e+01, 1.19449342e-12], [4.20000000e+01, 1.17400121e-12], [4.25000000e+01, 1.15413225e-12], [4.30000000e+01, 1.13478322e-12], [4.35000000e+01, 1.11598960e-12], [4.40000000e+01, 1.09764123e-12], [4.45000000e+01, 1.07962483e-12], [4.50000000e+01, 1.06202586e-12], [4.55000000e+01, 1.04493374e-12], [4.60000000e+01, 1.02821495e-12], [4.65000000e+01, 1.01187883e-12], [4.70000000e+01, 9.96036092e-13], [4.75000000e+01, 9.80360455e-13], [4.80000000e+01, 9.65162058e-13], [4.85000000e+01, 9.50366563e-13], [4.90000000e+01, 9.35695396e-13], [4.95000000e+01, 9.21522663e-13], [5.00000000e+01, 9.07535237e-13]])
    f = interpol.interp1d(Cm_mat[:,0], Cm_mat[:,1], "linear")
    return f(d)

def get_Cm_lin(d):
    # Dont use this function, this was just used for comparing linear and quadratic interpolation
    return np.vectorize(get_Cm_interpol_lin)(d)


def get_Lm_nonvec(d):

    Lm_mat = np.array([[1.00000000e+00, 1.50370318e-11], [2.00000000e+00, 1.28176495e-11], [3.00000000e+00, 1.12177246e-11], [4.00000000e+00, 9.96369927e-12], [5.00000000e+00, 8.94052811e-12], [6.00000000e+00, 8.08682762e-12], [7.00000000e+00, 7.36199016e-12], [8.00000000e+00, 6.73831905e-12], [9.00000000e+00, 6.19573028e-12], [1.00000000e+01, 5.72086080e-12], [1.10000000e+01, 5.30132999e-12], [1.20000000e+01, 4.92745225e-12], [1.30000000e+01, 4.59423231e-12], [1.40000000e+01, 4.29470251e-12], [1.50000000e+01, 4.02450889e-12], [1.60000000e+01, 3.77978074e-12], [1.70000000e+01, 3.55708594e-12], [1.80000000e+01, 3.35442782e-12], [1.90000000e+01, 3.16833752e-12], [2.00000000e+01, 2.99802731e-12]])

    if d<Lm_mat[0,0]:
        raise ValueError('No prediction for Lm, d is too small.')
    if d>=Lm_mat[-1,0]:
        raise ValueError('No prediction for Lm, d is too large.')

    # find between which two datapoints our requested d is
    index = 0
    while d>Lm_mat[index+1, 0]:
        index += 1

    # linear interpolation between points
    dx = Lm_mat[index+1, 0] - Lm_mat[index, 0]
    dydx = (Lm_mat[index+1, 1] - Lm_mat[index, 1])/dx
    
    Lm = (Lm_mat[index, 1]) + dydx*(d-Lm_mat[index,0])

    return Lm

def get_Lm_interpol(d):

    Lm_mat = np.array([[1.00000000e+00, 6.48279082e-08], [1.50000000e+00, 5.91620803e-08], [2.00000000e+00, 5.46059042e-08], [2.50000000e+00, 5.07751493e-08], [3.00000000e+00, 4.74689856e-08], [3.50000000e+00, 4.45670927e-08], [4.00000000e+00, 4.19868068e-08], [4.50000000e+00, 3.96693064e-08], [5.00000000e+00, 3.75741658e-08], [5.50000000e+00, 3.56681552e-08], [6.00000000e+00, 3.39235177e-08], [6.50000000e+00, 3.23217629e-08], [7.00000000e+00, 3.08434918e-08], [7.50000000e+00, 2.94744478e-08], [8.00000000e+00, 2.82055298e-08], [8.50000000e+00, 2.70238560e-08], [9.00000000e+00, 2.59208129e-08], [9.50000000e+00, 2.48892049e-08], [1.00000000e+01, 2.39235226e-08], [1.05000000e+01, 2.30158950e-08], [1.10000000e+01, 2.21623782e-08], [1.15000000e+01, 2.13593537e-08], [1.20000000e+01, 2.06000111e-08], [1.25000000e+01, 1.98844456e-08], [1.30000000e+01, 1.92062237e-08], [1.35000000e+01, 1.85652624e-08], [1.40000000e+01, 1.79552454e-08], [1.45000000e+01, 1.73779204e-08], [1.50000000e+01, 1.68288901e-08], [1.55000000e+01, 1.63063835e-08], [1.60000000e+01, 1.58083156e-08], [1.65000000e+01, 1.53340006e-08], [1.70000000e+01, 1.48809347e-08], [1.75000000e+01, 1.44484964e-08], [1.80000000e+01, 1.40364231e-08], [1.85000000e+01, 1.36400850e-08], [1.90000000e+01, 1.32622665e-08], [1.95000000e+01, 1.29004265e-08], [2.00000000e+01, 1.25532645e-08], [2.05000000e+01, 1.22198701e-08], [2.10000000e+01, 1.19006996e-08], [2.15000000e+01, 1.15937214e-08], [2.20000000e+01, 1.12984729e-08], [2.25000000e+01, 1.10152382e-08], [2.30000000e+01, 1.07420984e-08], [2.35000000e+01, 1.04803839e-08], [2.40000000e+01, 1.02267365e-08], [2.45000000e+01, 9.98387055e-09], [2.50000000e+01, 9.74802302e-09], [2.55000000e+01, 9.52226049e-09], [2.60000000e+01, 9.30321312e-09], [2.65000000e+01, 9.09252269e-09], [2.70000000e+01, 8.88823706e-09], [2.75000000e+01, 8.69204170e-09], [2.80000000e+01, 8.50115046e-09], [2.85000000e+01, 8.31718142e-09], [2.90000000e+01, 8.13885625e-09], [2.95000000e+01, 7.96673412e-09], [3.00000000e+01, 7.79952439e-09], [3.05000000e+01, 7.63797300e-09], [3.10000000e+01, 7.48138250e-09], [3.15000000e+01, 7.32920778e-09], [3.20000000e+01, 7.18201367e-09], [3.25000000e+01, 7.03953178e-09], [3.30000000e+01, 6.90099168e-09], [3.35000000e+01, 6.76644290e-09], [3.40000000e+01, 6.63576046e-09], [3.45000000e+01, 6.50923874e-09], [3.50000000e+01, 6.38628033e-09], [3.55000000e+01, 6.26648767e-09], [3.60000000e+01, 6.15031841e-09], [3.65000000e+01, 6.03744253e-09], [3.70000000e+01, 5.92764339e-09], [3.75000000e+01, 5.82027812e-09], [3.80000000e+01, 5.71628720e-09], [3.85000000e+01, 5.61497196e-09], [3.90000000e+01, 5.51651140e-09], [3.95000000e+01, 5.42011767e-09], [4.00000000e+01, 5.32692304e-09], [4.05000000e+01, 5.23564869e-09], [4.10000000e+01, 5.14666389e-09], [4.15000000e+01, 5.06060430e-09], [4.20000000e+01, 4.97615863e-09], [4.25000000e+01, 4.89370893e-09], [4.30000000e+01, 4.81333409e-09], [4.35000000e+01, 4.73507477e-09], [4.40000000e+01, 4.65890639e-09], [4.45000000e+01, 4.58409294e-09], [4.50000000e+01, 4.51142011e-09], [4.55000000e+01, 4.44036962e-09], [4.60000000e+01, 4.37109900e-09], [4.65000000e+01, 4.30324538e-09], [4.70000000e+01, 4.23723426e-09], [4.75000000e+01, 4.17229331e-09], [4.80000000e+01, 4.10897439e-09], [4.85000000e+01, 4.04726532e-09], [4.90000000e+01, 3.98681151e-09], [4.95000000e+01, 3.92774298e-09], [5.00000000e+01, 3.86983856e-09]])
    f = interpol.interp1d(Lm_mat[:,0]/1e6, Lm_mat[:,1], "quadratic")

    return f(d)

def get_Lm(d):
    return np.vectorize(get_Lm_interpol)(d)


# Demonstration
"""
Cm_mat = generate_Cm()
Lm_mat = generate_Lm()
d = np.linspace(1,49,10000)

#plt.plot(Cm_mat[:,0], Cm_mat[:,1]*1e12, marker=".", linestyle="none", label="simulated points")
#plt.plot(d, get_Cm(d/1e6)*1e12, label="quadratic")
#plt.plot(Lm_mat[:,0], Lm_mat[:,1]*1e9, marker=".", linestyle="none", label="simulated points")
#plt.plot(d, get_Lm(d/1e6)*1e9, label="quadratic")

plt.plot(d, get_Lm(d/1e6)/get_Cm(d/1e6), label="ratio")


#plt.plot(d, get_Cm_lin(d)*1e12, label="linear")
#plt.vlines(x=Cm_mat[:,0], ymin=-0.03, ymax=0, color="red", linewidth=0.5)
#plt.plot(d, get_Cm(d)*1e12-get_Cm_lin(d)*1e12, label="difference")
plt.legend()
plt.xlabel("d [um]")
#plt.ylabel("Mutual capacitance [pF]")
plt.show()
"""

