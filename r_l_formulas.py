import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


import common_formulas as cf

def calc_Rs(om, Ck, RL):
    return (1+(om**2)*(Ck**2)*(RL**2))/((om**2)*(Ck**2)*RL)

def calc_Cs(om, Ck, RL):
    return Ck/(1+om**2*Ck**2*RL**2)



