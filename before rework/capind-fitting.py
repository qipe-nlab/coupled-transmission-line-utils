import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


filename = "C:\\Users\\mluka\\work\\data\\capind\\capacitance_d_1-200um_l100um.csv"
data = np.genfromtxt(filename, delimiter=",", skip_header=5)
data=data[::2]
print(data)
plt.plot(data[:,0], data[:,2])
plt.show()

