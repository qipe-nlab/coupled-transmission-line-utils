import numpy as np
import scipy.integrate as integrate
from scipy.constants import epsilon_0

def mcpw_capacitance_matrix(d: list[float], eps_r: float):
    """Calculate capacitance matrix of multiconductor coplanar waveguide (mcpw)
    The formula is from  
    - G. Ghione "An efficient, CAD-oriented model for the characteristic parameters of multiconductor buses in high-speed digital GaAs ICs"  
      https://link.springer.com/article/10.1007/BF01673907

    Args:
        d (list[float]): width and gap of mcpw
                         GND-gap(d_0)-conductor(d_1)-...-conductor(d_{N-2})-gap(d_{N-1})-GND
        eps_r: relative permittivity
    """
    if len(d)%2 == 0:
        raise ValueError("d must be odd length")
        
    d = np.array([0]+list(d), dtype=float)
    t_boundary = d.cumsum() # even length
    N = len(t_boundary)//2 - 1
    
    def SC_denominator(t: float, t_b: np.ndarray):
        """Calculate the denominator of Schwarz-Christoffel mapping

        Args:
            t (float): t-plane coordinate
            t_b (np.ndarray): boundary coordinate between conductor and gap
        """
        D = 1
        for i in range(len(t_b)):
            D *= np.emath.sqrt(t-t_b[i])
        return complex(D)
    
    F = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            F[i][j] = integrate.quad(lambda x: x**i/SC_denominator(x, t_boundary), t_boundary[2*j+1], t_boundary[2*j+2], complex_func=True)[0]
    F = F.real
    
    G = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            G[i][j] = integrate.quad(lambda x: x**i/SC_denominator(x, t_boundary), t_boundary[2*j], t_boundary[2*j+1], complex_func=True)[0]
    G = G.imag
    G_inv = np.linalg.inv(G)
    
    U_inv = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            if i == j:
                U_inv[i][j] = 1
            elif i+1 == j:
                U_inv[i][j] = -1
    C = -epsilon_0*(1+eps_r)*U_inv@G_inv@F
    return C

# w = 5e-6
# s = 7.5e-6
# eps_r = 11.7

def calc_self_and_coupling_capacitance(d_list, w, s, eps_r):
    cs = []
    cm = []
    for d in d_list:
        p = [s, w, s, d, s, w, s]
        C = mcpw_capacitance_matrix(p, eps_r)
        cs.append(C[0][0])
        cm.append(C[0][2])
    return cs, cm