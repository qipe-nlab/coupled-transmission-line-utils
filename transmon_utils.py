## transmon functions
import numpy as np
from scipy import constants
from scipy.optimize import minimize, fsolve, root

import common_formulas as cf
# contants

e_charge = constants.e
hbar = constants.hbar
flux_q = hbar/(2*e_charge)

# functions

def EJ(L_J):

    val = flux_q**2/L_J

    return val

def EC(C_q):

    val = e_charge**2/(2*C_q)

    return val

def LJ(E_J):

    val = flux_q**2/E_J

    return val

def Cq(E_C):

    val = e_charge**2/(2*E_C)

    return val

def omega_q(EJ, EC):

    nu = (8*EJ/EC)**-0.5

    val = ((8*EJ*EC)**0.5 - EC*(1+nu+21*nu**2/2**3 + 608*nu**3/2**6 + 5219*nu**4/2**7))/hbar

    return val

def alpha(EJ, EC):

    nu = (8*EJ/EC)**-0.5

    val = -EC*(1+9*nu/2**2+81*nu**2/2**3 + 3645*nu**3/2**6 + 46899*nu**4/2**7)/hbar

    return val

def EC_from_spec(omega_q_val, EJ_val):

    if np.isnan(omega_q_val) or np.isnan(EJ_val):
        return np.nan

    fit_func = lambda Ec_fit: (omega_q(EJ_val, Ec_fit[0])-omega_q_val)**2

    Ec_initial = EJ_val/50
    
    #Note: for soome reason minimize does not yield a good solution
    res = fsolve(fit_func, Ec_initial)
    
    Ec_fitted = res
    return Ec_fitted[0]

def EJ_EC_from_spec(omega_q_val, alpha_val):

    ## returns values in MHz!

    if np.isnan(omega_q_val) or np.isnan(alpha_val):

        EJ_fitted, EC_fitted = (np.nan, np.nan)

    else:

        fit_func = lambda E_fit: (omega_q(E_fit[0]*hbar, E_fit[1]*hbar) - omega_q_val)**2 + (alpha(E_fit[0]*hbar, E_fit[1]*hbar) - alpha_val)**2

        E_initial = np.array([EJ(4e-9), EC(50e-15)])/hbar
        
        res = minimize(fit_func, E_initial, bounds = [(E_initial[0]*0.1,E_initial[0]*10), (E_initial[1]*0.1,E_initial[1]*10)])
        #res = root(fit_func, E_initial)
        EJ_fitted, EC_fitted = res.x

    return EJ_fitted*hbar, EC_fitted*hbar

def Cq_Lj_from_spec(omega_q_val, alpha_val):

    if np.isnan(omega_q_val) or np.isnan(alpha_val):

        Cq_fitted, Lj_fitted = (np.nan, np.nan)

    else:

        fit_func = lambda E_fit: (omega_q(E_fit[0]*hbar, E_fit[1]*hbar) - omega_q_val)**2 + (alpha(E_fit[0]*hbar, E_fit[1]*hbar) - alpha_val)**2 

        E_initial = np.array([EJ(4e-9), EC(50e-15)])/hbar

        res = minimize(fit_func, E_initial, bounds = [(E_initial[0]*0.1,E_initial[0]*10), (E_initial[1]*0.1,E_initial[1]*10)])

        EJ_fitted, EC_fitted = res.x

        Cq_fitted = Cq(EC_fitted*hbar)

        Lj_fitted = LJ(EJ_fitted*hbar)

    return Cq_fitted, Lj_fitted

def T1_radiative_spectrum(Cq_shunt, Z11):

    Y11 = 1/ Z11

    R_eff = 1/np.abs(np.real(Y11))

    eff_T1 = R_eff * Cq_shunt

    return eff_T1

## qubit-resonator coupling functions

def dispersive_shift(g, EC, omega_q, omega_r):

    ## function returns chi - i.e. half the full dispersive shift.
    # Expression is analytic approximation valid for Transmons

    delta = omega_q - omega_r

    val = -g**2*EC/hbar/(delta*(delta-EC/hbar))

    return val

def g(chi, EC, omega_q, omega_r):
    f = lambda x : dispersive_shift(x, EC, omega_q, omega_r) - chi

    chi_0 = -1 *2*np.pi*10**6

    val = np.abs(fsolve(f, chi_0))
    return val[0]

def g_from_r(r_coupling_param, omega_q, omega_r, C_shunt):

    # here r is a dimensional turns ratio element of a multiport Belevitch transformer.
    # the dimension is 1/(Farad)**0.5 
    ## see eq. 54 of https://arxiv.org/abs/1712.08154

    val = 0.5*r_coupling_param*(C_shunt*omega_q*omega_r)**0.5

    return val[0]

def r_from_g(g, omega_q, omega_r, C_shunt):

    f = lambda x : g_from_r(x, omega_q, omega_r, C_shunt) - g

    g_0 = 100*10**6*2*np.pi

    val = fsolve(f, g_0)

    return val[0]



fqs = [4.2e9, 4.7e9, 4.8e9, 4.3e9]
frs = [6.2e9, 6.52e9, 6.7e9, 6.36e9]



# spec

for i in range(4):
    print("----------")
    print("QUBIT ", i+1)
    f_q = fqs[i]
    om_q = f_q*2*np.pi

    f_r = frs[i]
    om_r = f_r*2*np.pi

    Lj_spec = 12.5e-9
    Ej_spec = EJ(Lj_spec)


    
    chi = -5e6*2*np.pi

    Ec_spec = EC_from_spec(om_q, Ej_spec)
    alpha_spec = alpha(Ej_spec, Ec_spec)
    Cq_result, Lj_result = Cq_Lj_from_spec(om_q, alpha_spec)

    
    print(f'alpha = {alpha_spec/1e6/2/np.pi:.4f} MHz')
    print(f'Cq = {Cq_result*1e15:.4f} fF')
    print(f'Ej/Ec = {Ej_spec/Ec_spec:.4f}')

    g_result = g(chi, Ec_spec, om_q, om_r)
    r_result = r_from_g(g_result, om_q, om_r, Cq_result)
    print(f'g/2pi = {(g_result/2/np.pi/1e6):.4f} MHz')
    print(f'r = {r_result/1e5:.4f} [1e5]')

    print(f"res_len = {cf.lambda_by_4_Ltot(f_r)*1e6:.1f} [um]")
