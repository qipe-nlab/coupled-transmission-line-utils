import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm
from scipy import constants
from scipy.optimize import minimize
import sys

def Zcap(Cval, omega):

    return -1j/(Cval*omega)

def Zind(Lval, omega):

    return 1j*(Lval*omega)

def Zopen(Z0, phi):

    return Z0/(1j*np.tan(phi))

def Zshort(Z0, phi):

    return 1j*Z0*np.tan(phi)

def Zinput(Z0, Zload, phi):

    return Z0*(Zload+1j*Z0*np.tan(phi))/(Z0+1j*Zload*np.tan(phi))

def Zpara(Z1, Z2):

    return 1/(1/Z1+1/Z2)

def Zres(Cr, Lr, omega, Rr = None):

    val = Zpara(Zcap(Cr, omega),Zind(Lr, omega))

    if Rr is not None:

        val = Zpara(val, Rr)

    return val

def Q_factor(Cr, Rr, omega):

    val = omega * Rr * Cr

    return val

def Zfilter(Cf, Cr, Lr, omega, R_purc = None):

    if R_purc is None:

        val = Zres(Cr,Lr, omega) + Zcap(Cf, omega)

    else:

        val = Zres(Cr,Lr, omega, R_purc) + Zcap(Cf, omega)

    return val

def C_eff(Z, omega):

    val = np.real(-1j/(Z*omega))

    return val

def res_omega(Cr, Lr):

    return 1/(Cr*Lr)**0.5

#####

## transmon functions

e_charge = constants.e
hbar = constants.hbar
flux_q = hbar/(2*e_charge)

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

def EJ_EC_from_spec(omega_q_val, alpha_val):

    if np.isnan(omega_q_val) or np.isnan(alpha_val):

        EJ_fitted, EC_fitted = (np.nan, np.nan)

    else:

        fit_func = lambda E_fit: (omega_q(E_fit[0]*hbar, E_fit[1]*hbar) - omega_q_val)**2 + (alpha(E_fit[0]*hbar, E_fit[1]*hbar) - alpha_val)**2 

        E_initial = np.array([EJ(4e-9), EC(50e-15)])/hbar

        res = minimize(fit_func, E_initial, bounds = [(E_initial[0]*0.5,E_initial[0]*2), (E_initial[1]*0.5,E_initial[1]*2)])

        EJ_fitted, EC_fitted = res.x 

    return EJ_fitted*hbar, EC_fitted*hbar

def Cq_Lj_from_spec(omega_q_val, alpha_val):

    if np.isnan(omega_q_val) or np.isnan(alpha_val):

        Cq_fitted, Lj_fitted = (np.nan, np.nan)

    else:

        fit_func = lambda E_fit: (omega_q(E_fit[0]*hbar, E_fit[1]*hbar) - omega_q_val)**2 + (alpha(E_fit[0]*hbar, E_fit[1]*hbar) - alpha_val)**2 

        E_initial = np.array([EJ(4e-9), EC(50e-15)])/hbar

        res = minimize(fit_func, E_initial, bounds = [(E_initial[0]*0.5,E_initial[0]*2), (E_initial[1]*0.5,E_initial[1]*2)])

        EJ_fitted, EC_fitted = res.x

        Cq_fitted = Cq(EC_fitted*hbar)

        Lj_fitted = LJ(EJ_fitted*hbar)

    return Cq_fitted, Lj_fitted

##### 

class circuit(object):

    def __init__(self, Zenv, Cbf1, Cbf2, Lbf, Cext, Cp, Lp, Cpr, Lpr, Cr, Lr, Cg, Cq, Lq, single_resonator = False):

        self.Zenv = Zenv
        self.Cbf1 = Cbf1
        self.Cbf2 = Cbf2
        self.Lbf = Lbf
        self.Cext = Cext
        self.Cp = Cp
        self.Lp = Lp
        self.Cpr = Cpr
        self.Lpr = Lpr
        self.Cr = Cr
        self.Lr = Lr
        self.Cg = Cg
        self.Cq = Cq
        self.Lq = Lq
        self.single_resonator = single_resonator

    def Z_qubit_port(self, omega):

        Z_qshunt = Zcap(self.Cq, omega) 

        Z_coupler_qr = Zcap(self.Cg, omega) 

        Z_r = Zres(self.Cr, self.Lr, omega)

        Z_coupler_pr = Zres(self.Cpr, self.Lpr, omega)

        Z_p = Zres(self.Cp, self.Lp, omega)

        Z_ext = Zcap(self.Cext, omega)

        Z_bf = Zfilter(self.Cbf1, self.Cbf2, self.Lbf, omega)

        if self.single_resonator:

            Zval = Zpara(Z_qshunt, Z_coupler_qr + Zpara(Z_p, Z_ext + Zpara(Z_bf, self.Zenv)))

        else:
            Zval = Zpara(Z_qshunt, Z_coupler_qr + Zpara(Z_r, Z_coupler_pr + Zpara(Z_p, Z_ext + Zpara(Z_bf, self.Zenv))))

        return Zval

class circuit_ham(object):

    def __init__(self, kbf, wbf1, wbf2, kp, wr, wrf, J, Zchar, g, wq, alpha_q, Zenv, single_resonator = False):

        self.kbf = kbf
        self.wbf1 = wbf1
        self.wbf2 = wbf2
        self.kp = kp
        self.wr = wr
        self.wrf = wrf
        self.J = J
        self.Zchar = Zchar
        self.g = g
        self.wq = wq
        self.alpha_q = alpha_q
        self.Zenv = Zenv
        self.single_resonator = single_resonator

    def get_qubit_circuit_params(self):

        ### test
        # Cq_val = 100e-15
        # Lq_val = 5e-9

        Cq_val, Lq_val = Cq_Lj_from_spec(wq, alpha_q)

        # print('Cq_val:', Cq_val)
        # print('Lq_val:', Lq_val)

        return Cq_val, Lq_val

    def get_bf_circuit_params(self):

        Cbf1_val = self.kbf/self.Zenv * 1/(self.wbf2)**2 * 1/((self.wbf2/self.wbf1)**2-1)

        Cbf2_val = self.kbf/self.Zenv * 1/(self.wbf2)**2 * 1/((self.wbf2/self.wbf1)**2-1)**2

        Lbf_val = self.Zenv/self.kbf * ((self.wbf2/self.wbf1)**2-1)**2 

        return Cbf1_val, Cbf2_val, Lbf_val

    def get_effective_readout_res_decay_rate(self):
        
        ### here valid for resonant double resonator        
        
        delta = 0

        kappa_eff_val = self.kp/2 * (1 - np.real(np.sqrt((self.kp - 2*1j*delta)**2 - 16*self.J**2))/self.kp)

        return kappa_eff_val

    def get_double_resonator_circuit_params(self):

        Cp_val = np.pi/(4*self.wr * self.Zchar)

        Lp_val = 4 * self.Zchar / (self.wr * np.pi)

        Cr_val = Cp_val
        Lr_val = Lp_val

        #### decay factor accounts for the shift to the Purcell resonator mode due to its external loss. This is necessary to keep the Purcell resonator and readout resonator on resonance

        numerator_factor = (self.kp/(4*Cp_val*self.Zenv))**0.5

        decay_factor = 1/(1 + numerator_factor/self.wr)

        Cp_val = Cp_val * decay_factor

        Lp_val = Lp_val * decay_factor

        #print('decay_factor:', decay_factor)

        #### g_factor accounts for the shift to the readout resonator due to its coupling to the qubit. This is also necessary to keep the Purcell resonator and readout resonator on resonance
        
        Cq_val, Lq_val = Cq_Lj_from_spec(wq, alpha_q)
        
        approx_Cg_val =  2 * (Cr_val * Cq_val / (self.wr * self.wq))**0.5 * self.g 

        if self.single_resonator:

            Cp_val = Cp_val  - (1/((1/Cq_val) + (1/approx_Cg_val)))
            Cr_val = Cp_val

            Cpr_val = np.nan

            Lpr_val = np.nan

        else:

            Cr_val = Cr_val  - (1/((1/Cq_val) + (1/approx_Cg_val)))

            Cpr_val = (Cp_val * Cr_val)**0.5/(self.wr - self.wrf) * self.J

            Lpr_val = 1/(self.wrf**2 * Cpr_val)

        return Cp_val, Lp_val, Cpr_val, Lpr_val, Cr_val, Lr_val

    def get_Cext(self):

        Cp_val, Lp_val, Cpr_val, Lpr_val, Cr_val, Lr_val = self.get_double_resonator_circuit_params()

        Zbf_val = Zfilter(*self.get_bf_circuit_params(), self.wr)

        A = (1/self.Zenv)**2 + (1/np.absolute(Zbf_val))**2

        #Cext_val = 1/ self.wr * 1/ (1/(1/(self.Zenv * A * Cp_val * self.kp) - 1/(self.Zenv * A)**2)**0.5 - 1/(A * np.absolute(Zbf_val)))  

        ### testing
        Cext_val = (Cp_val * self.kp * self.Zenv * A)**0.5 / self.wr

        return Cext_val

    def get_Cg(self):

        Cp_val, Lp_val, Cpr_val, Lpr_val, Cr_val, Lr_val = self.get_double_resonator_circuit_params()

        Cq_val, Lq_val = self.get_qubit_circuit_params()

        Cg_val =  2 * (Cr_val * Cq_val / (self.wr * self.wq))**0.5 * self.g 

        return Cg_val

    def circuit_params_from_ham(self):

        Zenv_val = self.Zenv

        single_resonator = self.single_resonator

        Cbf1_val, Cbf2_val, Lbf_val = self.get_bf_circuit_params()

        Cp_val, Lp_val, Cpr_val, Lpr_val, Cr_val, Lr_val = self.get_double_resonator_circuit_params()

        Cext_val = self.get_Cext()

        Cq_val, Lq_val = self.get_qubit_circuit_params()

        Cg_val = self.get_Cg()    

        return Zenv_val, Cbf1_val, Cbf2_val, Lbf_val, Cext_val, Cp_val, Lp_val, Cpr_val, Lpr_val, Cr_val, Lr_val, Cg_val, Cq_val, Lq_val, single_resonator

    def get_circuit_from_ham(self):

        #Zenv_val, Cbf1_val, Cbf2_val, Lbf_val, Cext_val, Cp_val, Lp_val, Cpr_val, Lpr_val, Cr_val, Lr_val, Cg_val, Cq_val, Lq_val = 

        print('params:', self.circuit_params_from_ham())

        circuit_class = circuit(*self.circuit_params_from_ham())

        return circuit_class

### Return circuit parameters from Hamiltonian parameters 

cmap = cm.get_cmap('Set2')

omega_vals = np.linspace(7, 11, 200) *1e9 *2*np.pi

kbf = 1e9 * 2*np.pi
wbf1 = 9e9 * 2 * np.pi
wbf2 = 10e9 * 2 * np.pi
kp = 100e6 * 2*np.pi
wr = 10.25e9 * 2*np.pi
wrf = 8.3e9 * 2*np.pi

J = 40e6 * 2* np.pi
Zchar = 50
g = 350e6 * 2 * np.pi
wq = 8.2e9 * 2* np.pi
alpha_q = -400e6 * 2*np.pi
Zenv = 50

circuit_ham_both_notch = circuit_ham(kbf, wbf1, wbf2, kp, wr, wrf, J, Zchar, g, wq, alpha_q, Zenv)

test_circuit_dist_notch = circuit_ham(kbf, wbf1, wbf2, kp, wr, 1e3*2*np.pi, J, Zchar, g, wq, alpha_q, Zenv)

test_circuit_no_notch = circuit_ham(kbf, 1e3*2*np.pi, wbf2, kp, wr, 1e3*2*np.pi, J, Zchar, g, wq, alpha_q, Zenv)

test_circuit_ham_single_resonator = circuit_ham(kbf, 1e3*2*np.pi, wbf2, kp, wr, 1e3*2*np.pi, J, Zchar, g, wq, alpha_q, Zenv, single_resonator = True)

#test_circuit_ham_single_resonator_no_backside_res = circuit_ham(kbf, 1e3*2*np.pi, 1e6*2*np.pi, kp, wr, 1e3*2*np.pi, J, Zchar, g, wq, alpha_q, Zenv, single_resonator = True)

test_circuit = circuit_ham_both_notch.get_circuit_from_ham()
test_circuit2 = test_circuit_dist_notch.get_circuit_from_ham()
test_circuit3 = test_circuit_no_notch.get_circuit_from_ham()
test_circuit4 = test_circuit_ham_single_resonator.get_circuit_from_ham()
#test_circuit5 = test_circuit_ham_single_resonator_no_backside_res.get_circuit_from_ham()

Z_readout = test_circuit.Z_qubit_port(omega_vals)
Y_readout = 1/Z_readout

Z_readout2 = test_circuit2.Z_qubit_port(omega_vals)
Y_readout2 = 1/Z_readout2

Z_readout3 = test_circuit3.Z_qubit_port(omega_vals)
Y_readout3 = 1/Z_readout3

Z_readout4 = test_circuit4.Z_qubit_port(omega_vals)
Y_readout4 = 1/Z_readout4

# Z_readout5 = test_circuit5.Z_qubit_port(omega_vals)
# Y_readout5 = 1/Z_readout5

print('Cq:', test_circuit.Cq)

plt.plot(omega_vals/(2*np.pi*1e9), (test_circuit.Cq/np.real(Y_readout)/(1e-6)), linestyle = '--', color = cmap(1), label = 'with backside & distributed filter')
plt.plot(omega_vals/(2*np.pi*1e9), (test_circuit.Cq/np.real(Y_readout2)/(1e-6)), linestyle = '--', color = cmap(2), label = 'with backside filter')
plt.plot(omega_vals/(2*np.pi*1e9), (test_circuit.Cq/np.real(Y_readout3)/(1e-6)), linestyle = '--', color = cmap(3), label = 'double resonator')
plt.plot(omega_vals/(2*np.pi*1e9), (test_circuit.Cq/np.real(Y_readout4)/(1e-6)), linestyle = '--', color = cmap(4), label = 'single resontor')
#plt.plot(omega_vals/(2*np.pi*1e9), (test_circuit.Cq/np.real(Y_readout5)/(1e-6)), linestyle = '--', color = cmap(5), label = 'single resontor without backside res')

plt.title('Radiative T1 limit for given mode frequencies and couplings')
plt.xlabel('qubit frequency (GHz)')
plt.ylabel('T1 limit (us)')
plt.yscale('log')
plt.legend(loc = 'lower left')
plt.show()
