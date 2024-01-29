import numpy as np
import matplotlib.pyplot as plt

def test_func(omega_n, omega_r, omega_p):

    fac_r = (omega_r/omega_n - omega_n/omega_r)
    fac_p = (omega_p/omega_n - omega_n/omega_p)

    val = 16/np.pi**2 * np.cos(np.pi * omega_n/ (2*omega_r)) * np.cos(np.pi * omega_n/ (2*omega_p)) / (fac_r*fac_p)

    return val

def test_func1(omega_r, omega):

    fac_r = (omega_r/omega - omega/omega_r)

    val = 4/np.pi * np.cos(np.pi * omega/ (2*omega_r)) / fac_r

    return val

omega_r = 10 # np.linspace( 4, 11, 51)
omega_rs = np.linspace( 4, 11, 51) #10
omega_p = 10.5

func_vals = test_func1(omega_r, omega_rs)

plt.plot(omega_rs, func_vals)
plt.show()