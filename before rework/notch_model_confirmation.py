import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from Coupled_transmission_lines_testing import find_notch_freq, find_lambda_freq


def get_error(l_Rf, l_Rn, l_Gf, l_Gn, l_c, plot=0):

    winding = np.pi*10e-6/4
    # d=20
    Lm = 451e-12/1e-3 #100e-12
    Cm = 3e-15/1e-3 #1e-15


    # read data
    sample_name = "coupled_" + str(int(l_Rf*1e6)) + "_" + str(int(l_Rn*1e6)) + "_" + str(int(l_Gf*1e6)) + "_" + str(int(l_Gn*1e6)) + "_" + str(int(l_c*1e6)) + "_20.csv"
    filename = "C:\\Users\\mluka\\work\\data\\coupled\\" + sample_name
    data = np.genfromtxt(filename,delimiter="," ,skip_header=5)

    # notches
    notch_freq = find_notch_freq(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm)
    # find minimum
    notch_index = np.argmin(data[:,2])

    #   l/2
    cpw_length1 = l_c + l_Gf + l_Gn + 4*winding
    cpw_length2 = l_c + l_Rf + l_Rn + 4*winding
    cpw1_freq = find_lambda_freq(cpw_length1)
    cpw2_freq = find_lambda_freq(cpw_length2)


    if plot!=0:
        plt.plot(data[:,0], np.log(data[:,2]))
        plt.plot(data[notch_index, 0], np.log(data[notch_index,2]), linestyle="none", marker="o", color="red")
        plt.vlines(x=notch_freq, ymax=5, ymin=-5, color="red")
        plt.vlines(x=[cpw1_freq, cpw2_freq], ymax=5, ymin=-5, color="green")
        plt.xlabel("f [GHz]")
        plt.ylabel("log(V2) [a.u.]")
        plt.show()

    
    return notch_freq-data[notch_index, 0]


l_Rf = 2e-3
l_Rn = 1.2e-3
l_Gf = 1.6e-3
l_Gn = 1.3e-3
l_c = 1.3e-3

# format, Rf, Rn, Gf, Gn, C
samples = [[2e-3, 1.2e-3, 1.6e-3, 1.3e-3, 1.3e-3], 
[1.4e-3, 0.9e-3, 1.6e-3, 1e-3, 1e-3], 
[1.4e-3, 1.2e-3, 1.6e-3, 1e-3, 0.7e-3], 
[1.7e-3, 0.9e-3, 1.3e-3, 1e-3, 1e-3],
[2e-3, 1.2e-3, 1e-3, 0.7e-3, 0.7e-3]]

errors = []

for i in range(len(samples)):
    errors.append(get_error(samples[i][0], samples[i][1], samples[i][2], samples[i][3], samples[i][4], plot=1))

print(errors)