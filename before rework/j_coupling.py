import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from Coupled_transmission_lines_testing import get_lumped_elements, Z_transfer_total


def par(Z1, Z2):
    return 1/(1/Z1 + 1/Z2)

def total_impedance(omega, C1, L1, C2, L2, Cg, Lg):
    Z1 = 1/(1j*omega*C1)
    Z2 = 1/(1j*omega*Cg + 1/(1j*omega*Lg))
    Z3 = 1/(1j*omega*C2)

    Z_tot = Z1 * Z3 / Z2 * 1/(1 + (Z1 + Z3)/Z2)
    return Z_tot

def j_coupling_formula(om1, om2, Z1, Z2, L1, L2):
    return -0.25*np.sqrt(om1*om2/L1/L2)*np.imag(Z1/om1+Z2/om2)

def j_coupling(om1, om2, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm):
    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm)
    Z1 = total_impedance(om1, C1, L1, C2, L2, Cg, Lg)
    Z2 = total_impedance(om2, C1, L1, C2, L2, Cg, Lg)
    j = j_coupling(om1, om2, Z1, Z2, L1, L2)
    return j

Lm = 451e-12/1e-3
Cm = 3e-15/1e-3
l_Rf = 1.7e-3
l_Rn = 0.9e-3
l_Gf = 1.3e-3
l_Gn = 1e-3
l_c = 1e-3


omegas = np.linspace(14*np.pi*1e9, 19*np.pi*1e9, 10000)

#Cm = 50e-12
#Lm = 65**2 * Cm * 0.45
def plot_model(Lm, Cm):
    

    C1, L1, C2, L2, Cg, Lg = get_lumped_elements(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm)
    print(C1, L1, C2, L2, Cg, Lg)
    Z1 = total_impedance(omegas, C1, L1, C2, L2, Cg, Lg)

    Zs = np.abs(Zs)
    Zs3 = Z_transfer_total(3*10**8/2.5, 65, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas)

    plt.plot(omegas/(2*np.pi*1e9), Zs, label="lumped element model  "+str(Lm)+"  "+str(Cm), color="blue")
    #plt.plot(omegas/(2*np.pi*1e9), np.log(Zs3), label="dist element model  "+str(Lm)+"  "+str(Cm), color="red")


filename = "C:\\Users\\mluka\\work\\data\\j_coupling_9.txt"
data = np.genfromtxt(filename, skip_header=5)
print(data)
difference = np.zeros(int(len(data)/2))
freqs = np.zeros(int(len(data)/2))
for i in range(len(difference)):
    difference[i] = np.abs(data[2*i, 1] - data[2*i+1, 1])
    freqs[i] = data[2*i, 0]

plt.plot(freqs, difference, linestyle="none", marker=".", color="red")

plt.xlabel("d [um]")
plt.ylabel("J-coupling")
plt.legend()
plt.show()

"""
filename = "C:\\Users\\mluka\\work\\data\\avoided_crossing_2.txt"
data = np.genfromtxt(filename, skip_header=5)
print(data)
plt.plot(data[:,0], data[:,1], linestyle="none", marker="o")
plt.show()
"""