
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as lin
import scipy.interpolate as interpol



######
# Cq #
######

def generate_Cq(filename):

    data = np.genfromtxt(filename,delimiter="," ,skip_header=5)
    len_data = int(np.shape(data)[0]/2)

    # Format: [[d,Cm],[d,Cm],[d,Cm]]
    Cq_mat = np.zeros((len_data, 2))

    for i in range(len_data):
        Cq_mat[i,0] = data[i*2,0]
        Cq_mat[i,1] = data[i*2,1]

    return Cq_mat

def get_Cq(d):
    f = interpol.interp1d(Cq_mat[:,0]/1e6, Cq_mat[:,1]/1e15, "quadratic", fill_value="extrapolate")
    return f(d)



######
# Cg #
######

def generate_Cg(filename):

    data = np.genfromtxt(filename,delimiter="," ,skip_header=5)
    len_data = int(np.shape(data)[0]/2)

    # Format: [[d,Cm],[d,Cm],[d,Cm]]
    Cg_mat = np.zeros((len_data, 2))

    for i in range(len_data):
        Cg_mat[i,0] = data[2*i,0]
        Cg_mat[i,1] = data[2*i,2]

    return Cg_mat

def get_Cg(d):
    f = interpol.interp1d(Cg_mat[:,0]/1e6, Cg_mat[:,1]/1e15, "linear", fill_value="extrapolate")
    return f(d)



######
# Cm #
######

def generate_Cm(filename):

    data = np.genfromtxt(filename,delimiter="," ,skip_header=5)
    len_data = int(np.shape(data)[0]/2)

    # Format: [[d,Cm],[d,Cm],[d,Cm]]
    Cm_mat = np.zeros((len_data, 2))

    for i in range(len_data):
        Cm_mat[i,0] = data[2*i,0]
        Cm_mat[i,1] = data[2*i,2]/1e-3 # account for simulation being of 1mm

    return Cm_mat


def get_Cm(d):
    f = interpol.interp1d(Cm_mat[:,0]/1e6, Cm_mat[:,1], "quadratic")
    return f(d)


######
# Lm #
######

def generate_Lm(filename):

    epsilon = 8.854e-12
    mu = 1.256e-6
    
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


def get_Lm(d):
    f = interpol.interp1d(Lm_mat[:,0]/1e6, Lm_mat[:,1], "quadratic")
    return f(d)

######
# Ck #
######

def generate_Ck(filename):

    data = np.genfromtxt(filename,delimiter="," ,skip_header=5)
    len_data = int(np.shape(data)[0]/2)

    # Format: [[d,Cm],[d,Cm],[d,Cm]]
    Ck_mat = np.zeros((len_data, 2))

    for i in range(len_data):
        Ck_mat[i,0] = data[2*i,0]
        Ck_mat[i,1] = data[2*i,2]

    return Ck_mat


def get_Ck(d):
    f = interpol.interp1d(Ck_mat[:,0], Ck_mat[:,1]/1e15, "linear", fill_value="extrapolate")
    return f(d)






##################
# Initialization #
##################


Cq_file = "C:\\Users\\mluka\\work\\dimensioning_algo\\cap_files\\cq_140-185.csv"
Cg_file = "C:\\Users\\mluka\\work\\dimensioning_algo\\cap_files\\cg_q-r-160_5-50.csv"

Cm_file = "C:\\Users\\mluka\\work\\dimensioning_algo\\cap_files\\capacitance_d-1-50-100_l-1000.csv"
Lm_file = "C:\\Users\\mluka\\work\\dimensioning_algo\\cap_files\\mcapacitance_d-1-50-100_l-1000.csv"

Ck_file = "C:\\Users\\mluka\\work\\dimensioning_algo\\cap_files\\ck_rd-20_5-75.csv"


Cg_mat = generate_Cg(Cg_file)
Cq_mat = generate_Cq(Cq_file)

Cm_mat = generate_Cm(Cm_file)
Lm_mat = generate_Lm(Lm_file)

Ck_mat = generate_Ck(Ck_file)





"""
d = np.linspace(1,150,10000)
plt.plot(Ck_mat[:,0], Ck_mat[:,1], marker=".", linestyle="none", label="simulated points")
plt.plot(d, get_Ck(d)*1e15, label="quadratic")
plt.legend()
plt.xlabel("via coupler deformation [um]")
plt.ylabel("Ck capacitance [fF]")
plt.show()
"""
"""

quit()
d = np.linspace(1,49,10000)

plt.plot(Cm_mat[:,0], Cm_mat[:,1]*1e12, marker=".", linestyle="none", label="simulated points")
plt.plot(d, get_Cm(d/1e6)*1e12, label="quadratic")
#plt.plot(Lm_mat[:,0], Lm_mat[:,1]*1e9, marker=".", linestyle="none", label="simulated points")
#plt.plot(d, get_Lm(d/1e6)*1e9, label="quadratic")

#plt.plot(d, get_Lm(d/1e6)/get_Cm(d/1e6), label="ratio")


#plt.plot(d, get_Cm_lin(d)*1e12, label="linear")
#plt.vlines(x=Cm_mat[:,0], ymin=-0.03, ymax=0, color="red", linewidth=0.5)
#plt.plot(d, get_Cm(d)*1e12-get_Cm_lin(d)*1e12, label="difference")
plt.legend()
plt.xlabel("d [um]")
#plt.ylabel("Mutual capacitance [pF]")
plt.show()



"""
"""
d = np.linspace(140e-6,200e-6,10000)
plt.plot(Cq_mat[:,0], Cq_mat[:,1], marker=".", linestyle="none", label="simulated points")
plt.plot(d*1e6, get_Cq(d)*1e15, label="quadratic")
plt.legend()
plt.xlabel("Qubit diameter [um]")
plt.ylabel("Qubit to ground capacitance [fF]")
plt.show()

"""
