# This file takes care of fitting the values obtained from the 2d simulations

######################
### IMPORTANT NOTE ###
######################

# The functions that should be used in other scripts are get_Cq and get_Lg

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as lin
import scipy.interpolate as interpol


def generate_Ck(sample_name="r_v_k_cap_10-70.csv"):

    filename = "C:\\Users\\mluka\\work\\data\\" + sample_name
    data = np.genfromtxt(filename,delimiter="," ,skip_header=5)
    len_data = int(np.shape(data)[0]/2)

    # Format: [[d,Cm],[d,Cm],[d,Cm]]
    Ck_mat = np.zeros((len_data, 2))

    for i in range(len_data):
        Ck_mat[i,0] = data[2*i,0]
        Ck_mat[i,1] = data[2*i,2]

    return Ck_mat

def get_Ck_interpol(d):
    #Ck_mat = np.array([[1., 19.19846918], [5.33333333, 19.06401538], [9.66666667, 19.17132499], [14., 19.26432611], [18.33333333, 19.26409409], [22.66666667, 19.22273698], [27., 19.10885712], [31.33333333, 18.87609771], [35.66666667, 18.64365432], [40., 18.4633244]])
    Ck_mat = np.array([[10., 8.52991762], [16.66666667, 10.20708956], [23.33333333, 11.95607729], [30., 13.70097146], [36.66666667, 15.50732765], [43.33333333, 17.38219853], [50., 19.24260636], [56.66666667, 21.03209343], [63.33333333, 22.92562689], [70., 24.84108554]])
    
    f = interpol.interp1d(Ck_mat[:,0], Ck_mat[:,1]/1e15, "quadratic")
    return f(d)


def get_Ck(d):
    return np.vectorize(get_Ck_interpol)(d)


#print(generate_Ck())
#quit()

#Ck_mat = generate_Ck()
#d = np.linspace(10,70,10000)
"""
plt.plot(Ck_mat[:,0], Ck_mat[:,1], marker=".", linestyle="none", label="simulated points")
plt.plot(d, get_Ck(d)*1e15, label="quadratic")
plt.legend()
plt.xlabel("via coupler deformation [um]")
plt.ylabel("Ck capacitance [fF]")
plt.show()
"""
#print(get_Ck(60))
