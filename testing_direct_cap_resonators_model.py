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

### direct lambda/4 resonators model

run_model_A = True
run_model_B = False
run_model_C = False

if run_model_A:

    phase_vel = 119919602
    Z0 = 65 #65

    Cl = 1/(default_phase_vel * default_Z0)
    Ll = default_Z0/default_phase_vel

    Cm = 5e-15

    l_Rf = 1.85e-3 #-0.51e-3
    l_Rn = 1.05e-3 #-0.51e-3
    l_Gf = 2.5e-3 #-0.51e-3
    l_Gn = 0.21e-3 #- 0.51e-3

    omegas = np.linspace(1, 20, 50) * 2*np.pi * 1e9

    Z_transfer_exact = Z_transfer_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65)
    C1, L1, C2, L2, Cg, Lg = get_lumped_elements_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65)
    Z_transfer_lumped_element = lumped_model_Z_transmission(omegas, C1, L1, C2, L2, Cg, Lg)

    plt.plot(omegas/(2*np.pi*1e9), np.abs(np.imag(Z_transfer_exact)), color = 'r')
    plt.plot(omegas/(2*np.pi*1e9), np.abs(np.imag(Z_transfer_lumped_element)), color = 'g')
    plt.yscale('log')
    plt.show()

    Cm_vals = np.linspace(0.25, 5, 40)*1e-15

    J_vals_direct_cap = np.array([J_coupling_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, Cm_v, phase_vel=phase_vel, Z0=Z0) for Cm_v in Cm_vals])
    J_vals_direct_cap_analytic = np.array([J_coupling_direct_cap_analytic(l_Gf, l_Gn, l_Rf, l_Rn, Cm_v, phase_vel=phase_vel, Z0=Z0) for Cm_v in Cm_vals])

    plt.plot(Cm_vals*1e15, J_vals_direct_cap/(2*np.pi*1e6), color = 'r')
    plt.plot(Cm_vals*1e15, J_vals_direct_cap_analytic/(2*np.pi*1e6), color = 'g')
    plt.yscale('log')
    plt.show()

    print('J_vals_direct_cap:', J_vals_direct_cap)
    print('J_vals_direct_cap_analytic:', J_vals_direct_cap_analytic)

### direct lambda/4 to lambda/2 resonators model

if run_model_B:

    phase_vel = 119919602
    Z0 = 65 #65

    Cl = 1/(default_phase_vel * default_Z0)
    Ll = default_Z0/default_phase_vel

    Cm = 5e-15

    l_Rf = 2.2e-3*2 #-0.51e-3
    l_Rn = 0.5e-3*2 #-0.51e-3
    l_Gf = 2.2e-3 #-0.51e-3
    l_Gn = 0.5e-3 #- 0.51e-3

    tau_Rf = l_Rf/phase_vel

    notch_freq_symbolic = np.pi/(2*tau_Rf)

    omegas = np.linspace(1, 20, 400) * 2*np.pi * 1e9

    Z_transfer_exact = Z_transfer_direct_cap_exact_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65)
    C1, L1, C2, L2, Cg, Lg = get_lumped_elements_direct_cap_B(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65)
    Z_transfer_lumped_element = lumped_model_Z_transmission(omegas, C1, L1, C2, L2, Cg, Lg)

    plt.plot(omegas/(2*np.pi*1e9), np.abs(np.imag(Z_transfer_exact)), color = 'r')
    plt.plot(omegas/(2*np.pi*1e9), np.abs(np.imag(Z_transfer_lumped_element)), color = 'g')
    plt.vlines(notch_freq_symbolic/(2*np.pi*1e9), 0.01, 10, color = 'b', linestyle = '--')

    plt.yscale('log')
    plt.show()

    CJs = np.linspace(0.2, 5, 10) * 1e-15

    l_Rf = 1.36e-3 * 2
    l_Rn = 1.35e-3 * 2

    J_vals_direct_cap = np.array([J_coupling_direct_cap_B(l_Gf, l_Gn, l_Rf, l_Rn, CJ, phase_vel=phase_vel, Z0=Z0) for CJ in CJs])
    J_vals_direct_symbolic = np.array([J_coupling_direct_cap_analytic_B(l_Gf, l_Gn, l_Rf, l_Rn, CJ, phase_vel=phase_vel, Z0=Z0) for CJ in CJs])

    print('J_vals_direct_cap:', J_vals_direct_cap)
    print('J_vals_direct_symbolic:', J_vals_direct_symbolic)

    plt.plot(CJs * 1e15, J_vals_direct_cap/(2*np.pi*1e6), color = 'r')
    plt.plot(CJs * 1e15, J_vals_direct_symbolic/(2*np.pi*1e6), color = 'g')
    plt.show()

    l_Gn_vals = np.linspace(0.1, 2.6, 50)*1e-3

    CJ = 2e-15

    J_vals_direct_cap = np.array([J_coupling_direct_cap_B(2.7e-3 - l_Gn_val, l_Gn_val, l_Rf, l_Rn, CJ, phase_vel=phase_vel, Z0=Z0) for l_Gn_val in l_Gn_vals])
    J_vals_direct_symbolic = np.array([J_coupling_direct_cap_analytic_B(2.7e-3 - l_Gn_val, l_Gn_val, l_Rf, l_Rn, CJ, phase_vel=phase_vel, Z0=Z0) for l_Gn_val in l_Gn_vals])

    print('J_vals_direct_cap:', J_vals_direct_cap)
    print('J_vals_direct_symbolic:', J_vals_direct_symbolic)

    plt.plot(l_Gn_vals * 1e3, J_vals_direct_cap/(2*np.pi*1e6), color = 'r')
    plt.plot(l_Gn_vals * 1e3, J_vals_direct_symbolic/(2*np.pi*1e6), color = 'g')
    plt.show()

    l_Rn_vals = np.linspace(0.125, 1.5, 15)*1e-3 * 2

    CJ = 2e-15

    J_vals_direct_cap = np.array([J_coupling_direct_cap_B(l_Gf, l_Gn, 5.2*1e-3 - l_Rn_val, l_Rn_val, CJ, phase_vel=phase_vel, Z0=Z0) for l_Rn_val in l_Rn_vals])
    J_vals_direct_symbolic = np.array([J_coupling_direct_cap_analytic_B(l_Gf, l_Gn, 5.2*1e-3 - l_Rn_val, l_Rn_val, CJ, phase_vel=phase_vel, Z0=Z0) for l_Rn_val in l_Rn_vals])

    print('J_vals_direct_cap:', J_vals_direct_cap)
    print('J_vals_direct_symbolic:', J_vals_direct_symbolic)

    plt.plot(l_Rn_vals * 1e3, J_vals_direct_cap/(2*np.pi*1e6), color = 'r')
    plt.plot(l_Rn_vals * 1e3, J_vals_direct_symbolic/(2*np.pi*1e6), color = 'g')
    plt.show()

### direct lambda/2 to lambda/2 resonators model

if run_model_C:

    phase_vel = 119919602
    Z0 = 65 #65

    Cl = 1/(default_phase_vel * default_Z0)
    Ll = default_Z0/default_phase_vel

    Cm = 5e-15

    l_Rf = 1.85e-3*2 #-0.51e-3
    l_Rn = 0.5e-3*2 #-0.51e-3
    l_Gf = 2.2e-3*2 #-0.51e-3
    l_Gn = 0.51e-3*2 #- 0.51e-3

    omegas = np.linspace(1, 20, 100) * 2*np.pi * 1e9

    Z_transfer_exact = Z_transfer_direct_cap_exact_C(l_Gf, l_Gn, l_Rf, l_Rn, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65)
    C1, L1, C2, L2, Cg1, Lg1, Cg2, Lg2 = get_lumped_elements_direct_cap_C(l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65)
    Z_transfer_lumped_element = lumped_model_Z_transmission(omegas, C1, L1, C2, L2, Cg1, Lg1, Cg2, Lg2)

    plt.plot(omegas/(2*np.pi*1e9), np.abs(np.imag(Z_transfer_exact)), color = 'r')
    plt.plot(omegas/(2*np.pi*1e9), np.abs(np.imag(Z_transfer_lumped_element)), color = 'g')
    plt.yscale('log')
    plt.show()