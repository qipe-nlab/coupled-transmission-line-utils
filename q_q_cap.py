# This file takes care of fitting the values obtained from the 2d simulations

######################
### IMPORTANT NOTE ###
######################

# The functions that should be used in other scripts are get_Cq and get_Cg

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as lin
import scipy.interpolate as interpol


def generate_Cc(sample_name="q_q_c_cap_3-45.csv"):

    filename = "C:\\Users\\mluka\\work\\data\\" + sample_name
    data = np.genfromtxt(filename,delimiter="," ,skip_header=5)
    len_data = int(np.shape(data)[0]/2)

    # Format: [[d,Cm],[d,Cm],[d,Cm]]
    Cc_mat = np.zeros((len_data, 2))

    for i in range(len_data):
        Cc_mat[i,0] = data[2*i,0]
        Cc_mat[i,1] = data[2*i,2]

    return Cc_mat

def get_Cc_interpol(angle):
    Cc_mat = np.array([[ 3., 3.78191815], [ 6., 3.81092934], [ 9., 4.00555599], [12., 4.2024138 ], [15., 4.39445769], [18., 4.59752178], [21., 4.78704021], [24., 4.98874791], [27., 5.19075928], [30., 5.38697242], [33., 5.60000171], [36., 5.79485923], [39., 5.99967028], [42., 6.20611702], [45., 6.41441731]])
    
    f = interpol.interp1d(Cc_mat[:,0], Cc_mat[:,1]/1e15, "quadratic")
    return f(angle)

def get_Cc_2dinterpol(qubit_d, angle):
    Ds = np.array([140, 145, 160])
    angles = np.array([3, 17, 31, 45])
    Cc = np.array([[3.82, 3.83, 3.833],[4.53, 4.58, 4.75],[5.43, 5.54, 5.85],[6.37, 6.52, 6.97]])

    f = interpol.interp2d(Ds,angles,Cc, bounds_error=True)
    return f(qubit_d*1e6, angle)/1e15
    
def get_Cc(qubit_d, angle):
    return np.vectorize(get_Cc_2dinterpol)(qubit_d, angle)



# Demonstration

#print(generate_Cc())

"""
Cc_mat = generate_Cc()
angle = np.linspace(3,45,1000)

Cc_mat = np.array([[3.74, 3.78, 3.80],[4.32, 4.52, 4.7],[5.1, 5.46, 5.78],[5.9, 6.41, 6.88]])
angles = np.array([3, 17, 31, 45])

plt.plot(angle, get_Cc(130e-6, angle)*1e15, label="130")
plt.plot(angle, get_Cc(145e-6, angle)*1e15, label="145")
plt.plot(angle, get_Cc(160e-6, angle)*1e15, label="160")
plt.plot(angle, get_Cc(150e-6, angle)*1e15, label="150")
plt.plot(angle, get_Cc(140e-6, angle)*1e15, label="140")

plt.plot(angles, Cc_mat[:,0], marker="o", linestyle="none", label="simulated points", color="red")
plt.plot(angles, Cc_mat[:,1], marker="o", linestyle="none", color="red")
plt.plot(angles, Cc_mat[:,2], marker="o", linestyle="none", color="red")

plt.legend()
plt.xlabel("Qubit diameter [um]")
plt.ylabel("Qubit to ground capacitance [fF]")
plt.show()
"""
