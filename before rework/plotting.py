import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


filename = "C:\\Users\\mluka\\work\\data\\j_coupling_10.txt"
data = np.genfromtxt(filename, skip_header=5)
print(data)
difference = np.zeros(int(len(data)/2))
freqs = np.zeros(int(len(data)/2))
for i in range(len(difference)):
    difference[i] = np.abs(data[2*i, 1] - data[2*i+1, 1])
    freqs[i] = data[2*i, 0]

plt.plot(freqs, difference, linestyle="none", marker="o", label="rougher mesh")

filename = "C:\\Users\\mluka\\work\\data\\j_coupling_8.txt"
data = np.genfromtxt(filename, skip_header=5)
print(data)
difference = np.zeros(int(len(data)/2))
freqs = np.zeros(int(len(data)/2))
for i in range(len(difference)):
    difference[i] = np.abs(data[2*i, 1] - data[2*i+1, 1])
    freqs[i] = data[2*i, 0]

#plt.plot(freqs, difference, linestyle="none", marker="o", color="green")

filename = "C:\\Users\\mluka\\work\\data\\j_coupling_9.txt"
data = np.genfromtxt(filename, skip_header=5)
print(data)
difference = np.zeros(int(len(data)/2))
freqs = np.zeros(int(len(data)/2))
for i in range(len(difference)):
    difference[i] = np.abs(data[2*i, 1] - data[2*i+1, 1])
    freqs[i] = data[2*i, 0]

plt.plot(freqs, difference, linestyle="none", marker=".", color="red")

plt.xlabel("f [GHz]")
plt.ylabel("Z21 [$ \Omega $]")
plt.yscale('log')
plt.legend()
plt.show()

"""
filename = "C:\\Users\\mluka\\work\\data\\avoided_crossing_2.txt"
data = np.genfromtxt(filename, skip_header=5)
print(data)
plt.plot(data[:,0], data[:,1], linestyle="none", marker="o")
plt.show()
"""
