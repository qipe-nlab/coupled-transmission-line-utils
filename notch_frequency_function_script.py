### example solutions

from exact_coupled_transmission_line_eqn_solver import *

### param set 1: for high freq range

phase_vel = 3*10**8/2.5
Z0 = 65

# Lm = 1.31e-9*1e1*1.8
# Cm = 5e-15*1e3*1.8
# l_Rf = 1.2e-3
# l_Rn = 1.2e-3
# l_Gf = 1.2e-3
# l_Gn = 1.2e-3
# l_c = 0.5e-3

## param high v2

Lm = 1.31e-9*1e1*2.2 
Cm = 5e-15*1e3*2.2
l_Rf = 0.99e-3 
l_Rn = 1.5e-3
l_Gf = 1.54e-3
l_Gn = 0.95e-3 
l_c = 0.5e-3

omega = 8 * 2*np.pi*1e9
C = 100 * 1e-15

L = 1/(C*omega**2)

print('L:', L)

L = 10 * 1e-9

C = 1/(L*omega**2)

print('C:', C)

### param set 2: for low freq range

# scale_fac = 1.6

# Lm = 1.31e-9*1e1*2.2 
# Cm = 5e-15*1e3*2.2
# l_Rf = 1.15e-3 * scale_fac
# l_Rn = 1.45e-3 * scale_fac
# l_Gf = 1.7e-3 * scale_fac
# l_Gn = 0.9e-3 * scale_fac
# l_c = 0.5e-3 * scale_fac

print('Lm:', Lm)
print('Cm:', Cm)
print('l_Rf:', l_Rf)
print('l_Rn:', l_Rn)
print('l_Gf:', l_Gf)
print('l_Gn:', l_Gn)
print('l_c:', l_c)
#sys.exit()

l_delt = 0.03e-3

l_Gn = l_Gn - l_delt
l_Gf = l_Gf + l_delt
l_Rn = l_Rn - l_delt
l_Rf = l_Rf + l_delt

test_notch_freq = find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=3*10**8/2.5, Z0=65)

print('notch_freq_exact (GHz):', test_notch_freq / (2*np.pi*1e9))

test_notch_freq_analytic = find_notch_filter_frequency_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=3*10**8/2.5, Z0=65, min_search=3*2*np.pi*10**9, max_search=14*2*np.pi*10**9, search_spacing=(0.5*2*np.pi*10**6))

print('notch_freq_analytic (GHz):', test_notch_freq_analytic / (2*np.pi*1e9))

notch_freq_agreement_percent = 100*(test_notch_freq_analytic - test_notch_freq )/ test_notch_freq

print('notch_freq_agreement_percent:', notch_freq_agreement_percent)
