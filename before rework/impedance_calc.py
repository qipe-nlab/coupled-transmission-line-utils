import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from coupled_transmission_line_util import get_lumped_elements, Z_transfer_total, lumped_element_Z21
from mutual_LC_util import get_Lm, get_Cm


def par(Z1, Z2):
    return 1/(1/Z1 + 1/Z2)

def total_impedance(omega, C1, L1, C2, L2, Cg, Lg):
    Z1 = 1/(1j*omega*C1 + 1/(1j*omega*L1))
    Z2 = 1/(1j*omega*Cg + 1/(1j*omega*Lg))
    Z3 = 1/(1j*omega*C2 + 1/(1j*omega*L2))

    Z_tot = Z1 * Z3 / Z2 * 1/(1 + (Z1 + Z3)/Z2)
    return Z_tot



#Lm = 1.25e-8
#Cm = 3e-12
#Cm = 2.9963e-15/1e-3
#Lm = 2.00692348e-09

l_Gn = 804e-6
l_Gf = 1653e-6
l_c = 374e-6
l_Rn = 827e-6
l_Rf = 1297e-6

d = 3.5 #in um

target = [9, 10, 8]

Lm = get_Lm(d)
Cm = get_Cm(d)
omegas = np.linspace(12*np.pi*1e9, 22*np.pi*1e9, 10000)

def plot_model(Lm, Cm):
    calibration = 600e-6
    #C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=119919602, Z0=64.6)
    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf+calibration/2, l_Gn+calibration/2, l_Rf+calibration/2, l_Rn+calibration/2, Lm, Cm)
    Zs = lumped_element_Z21(omegas, C1, L1, C2, L2, Cg, Lg)
    Zs = np.abs(Zs)

    plt.plot(omegas/(2*np.pi*1e9), Zs, label="lumped element model", color="blue")
    #plt.plot(omegas/(2*np.pi*1e9), np.log(Zs3), label="dist element model  "+str(Lm)+"  "+str(Cm), color="red")



#sample_name = "coupled_" + str(int(l_Rf*1e6)) + "_" + str(int(l_Rn*1e6)) + "_" + str(int(l_Gf*1e6)) + "_" + str(int(l_Gn*1e6)) + "_" + str(int(l_c*1e6)) + "_" + str(int(d)) + ".csv"

sample_name = "real_model_" + str(int(l_Gn*1e6)) + "_" + str(int(l_Gf*1e6)) + "_" + str(int(l_c*1e6)) + "_" + str(int(l_Rn*1e6)) + "_" + str(int(l_Rf*1e6)) + "_" + str(int(d)) + ".csv"

filename = "C:\\Users\\mluka\\work\\data\\coupled\\" + sample_name
data = np.genfromtxt(filename,delimiter="," ,skip_header=5)

#sample_name = "coupled_" + str(int(l_Rf*1e6)) + "_" + str(int(l_Rn*1e6)) + "_" + str(int(l_Gf*1e6)) + "_" + str(int(l_Gn*1e6)) + "_" + str(int(l_c*1e6)) + "_20_airbridge.csv"
#filename = "C:\\Users\\mluka\\work\\data\\coupled\\" + sample_name
#data2 = np.genfromtxt(filename,delimiter="," ,skip_header=5)


#sample_name = "coupled_" + str(int(l_Rf*1e6)) + "_" + str(int(l_Rn*1e6)) + "_" + str(int(l_Gf*1e6)) + "_" + str(int(l_Gn*1e6)) + "_" + str(int(l_c*1e6)) + "_20_via_ef.csv"
#filename = "C:\\Users\\mluka\\work\\data\\coupled\\" + sample_name
#data3 = np.genfromtxt(filename,delimiter="," ,skip_header=5)

plot_model(Lm, Cm)
#plt.plot(data[:,0], data[:,3], label="simulation", color="blue", marker=".")
#plt.vlines([target[0],target[1]], np.min(data[:,3])/100, np.max(data[:,3])*100, label="targeted resonators", color="orange")
#plt.vlines([target[2]], np.min(data[:,3])/100, np.max(data[:,3])*100, label="targeted notch", color="red")

plt.legend()
plt.xlabel("f [GHz]")
plt.ylabel("Z21 [$ \Omega $]")
plt.yscale('log')
plt.show()
