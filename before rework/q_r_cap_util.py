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


def generate_Cq(sample_name="q_r_gnd_cap_140-160.csv"):

    filename = "C:\\Users\\mluka\\work\\data\\" + sample_name
    data = np.genfromtxt(filename,delimiter="," ,skip_header=5)
    len_data = int(np.shape(data)[0]/2)

    # Format: [[d,Cm],[d,Cm],[d,Cm]]
    Cq_mat = np.zeros((len_data, 2))

    for i in range(len_data):
        Cq_mat[i,0] = data[2*i,0]
        Cq_mat[i,1] = data[2*i,1]

    return Cq_mat

def get_Cq_interpol(d):
    Cq_mat = np.array([[140.,37.50896241], [141.,37.9485321], [142.,38.38854591], [143.,38.83138371], [144.,39.26659103], [145.,39.70977016], [146.,40.17477219], [147.,40.63434331], [148.,41.10883209], [149.,41.58824328], [150.,42.07266291], [151.,42.55748014], [152.,43.04796863], [153.,43.56412635], [154.,44.08454274], [155.,44.62300469], [156.,45.13248794], [157.,45.68158668], [158.,46.23124304], [159.,46.78297593], [160.,47.35421507], [161.,47.97108699], [162.,48.54500841], [163.,49.16081849], [164.,49.7811336 ], [165.,50.41631571], [166.,51.06214908], [167.,51.74493168], [168.,52.45773182], [169.,53.14300458], [170.,53.895217]])

    
    f = interpol.interp1d(Cq_mat[:,0]/1e6, Cq_mat[:,1]/1e15, "quadratic")
    return f(d)


def get_Cq(d):
    return np.vectorize(get_Cq_interpol)(d)


def generate_Cg(sample_name="q_r_g_cap_10-60.csv"):

    filename = "C:\\Users\\mluka\\work\\data\\" + sample_name
    data = np.genfromtxt(filename,delimiter="," ,skip_header=5)
    len_data = int(np.shape(data)[0]/2)

    # Format: [[d,Cm],[d,Cm],[d,Cm]]
    Cg_mat = np.zeros((len_data, 2))

    for i in range(len_data):
        Cg_mat[i,0] = data[2*i,0]
        Cg_mat[i,1] = data[2*i,2]

    return Cg_mat

def get_Cg_interpol(d):
    Cg_mat = np.array([[10.0, 5.5108818], [11.57894737, 5.6540051], [13.15789474, 5.80858648], [14.73684211, 5.96304858], [16.31578947, 6.11525188], [17.89473684, 6.27937674], [19.47368421, 6.43997494], [21.05263158, 6.59937484], [22.63157895, 6.75665746], [24.21052632, 6.90972162], [25.78947368, 7.06415122], [27.36842105, 7.21486288], [28.94736842, 7.3637014], [30.52631579, 7.50938626], [32.10526316, 7.65826108], [33.68421053, 7.75570041], [35.26315789, 7.86180543], [36.84210526, 7.96536446], [38.42105263, 8.06535024], [40.0, 8.13610712]])

    #Cg_mat = np.array([[10., 5.98485942], [15.55555556, 6.52666317], [21.11111111, 7.08532197], [26.66666667, 7.6379622 ], [32.22222222, 8.14729809], [37.77777778, 8.55324916], [43.33333333, 8.8418195 ], [48.88888889, 9.06617667], [54.44444444, 9.24371826], [60.0, 9.43441922]])
    
    f = interpol.interp1d(Cg_mat[:,0]/1e6, Cg_mat[:,1]/1e15, "quadratic")
    return f(d)


def get_Cg(d):
    return np.vectorize(get_Cg_interpol)(d)


"""
print(generate_Cg())


Cg_mat = generate_Cg()
d = np.linspace(10e-6,59e-6,10000)

plt.plot(Cg_mat[:,0], Cg_mat[:,1], marker=".", linestyle="none", label="simulated points")
plt.plot(d*1e6, get_Cg(d)*1e15, label="quadratic")
plt.legend()
plt.xlabel("Qubit diameter [um]")
plt.ylabel("Qubit to ground capacitance [fF]")
plt.show()
"""
