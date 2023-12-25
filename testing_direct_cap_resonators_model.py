import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import sys
import cap_util as cap

default_phase_vel = 119919602
default_Z0 = 65

from exact_coupled_transmission_line_eqn_solver import *

# omegas = np.linspace(2, 12, 100) * 2*np.pi * 1e9

#phase_vel=3*10**8/2.5
#Z0 = 65

phase_vel = 119919602
Z0 = 65 #65

Cl = 1/(default_phase_vel * default_Z0)
Ll = default_Z0/default_phase_vel

Cm = 3e-15

l_Rf = 1.85e-3*2 #-0.51e-3
l_Rn = 1.05e-3*2 #-0.51e-3
l_Gf = 2.2e-3 #-0.51e-3
l_Gn = 0.51e-3 #- 0.51e-3

omegas = np.linspace(1, 20, 50) * 2*np.pi * 1e9

Z_transfer_exact = Z_transfer_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65)
C1, L1, C2, L2, Cg, Lg = get_lumped_elements_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65)
Z_transfer_lumped_element = lumped_model_Z_transmission(omegas, C1, L1, C2, L2, Cg, Lg)

plt.plot(omegas/(2*np.pi*1e9), np.abs(np.imag(Z_transfer_exact)), color = 'r')
plt.plot(omegas/(2*np.pi*1e9), np.abs(np.imag(Z_transfer_lumped_element)), color = 'g')
plt.yscale('log')
plt.show()

###
