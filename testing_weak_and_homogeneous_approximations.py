import seaborn as sns
import numpy as np
import scipy as sp
from scipy import constants
import matplotlib.pyplot as plt

from matplotlib.colors import ListedColormap

import sys
import cap_util as cap

hex_codes = sns.color_palette().as_hex()
hex_codes2 = sns.color_palette('husl', 9).as_hex()
my_cmap = ListedColormap(hex_codes)
my_cmap2 = ListedColormap(hex_codes2)
my_cmap3 = ListedColormap(["#9b59b6", "#34495e"])

d_vals = np.linspace(1, 40, 20) * 1e-6

# print('d_vals:', d_vals)

Cm_vals = cap.get_Cm(d_vals)
Lm_vals = cap.get_Lm(d_vals)
Zm_vals = cap.get_Zm(d_vals)

################

default_phase_vel = constants.c / np.sqrt(6.351)
default_Z0 = 65.62

Cl = 1/(default_phase_vel * default_Z0)
Ll = default_Z0/default_phase_vel

################

weak_coupling_validity = Lm_vals / Ll

plt.plot(d_vals* 1e6, weak_coupling_validity, color = my_cmap3(1), marker = 'o', markersize = 8, linestyle = 'None')

ax = plt.gca()

# Customize the border settings
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2.5)
ax.spines['left'].set_linewidth(2.5)

# Adjust tick thickness and label size
ax.tick_params(axis='both', which='both', width=2.5)
ax.tick_params(axis='both', labelsize=25)

ax.set_xticks([0, 5, 10, 15, 20, 25, 30, 35, 40 ])
ax.set_xticklabels(['0', '5', '10', '15', '20', '25', '30', '35', '40'])

# Set axis labels and font size
ax.set_xlabel(r'$d$ (um)', fontsize=30)
ax.set_ylabel(r'$l_m/l_c$', fontsize=30)

# Add legend with fontsize 20

#ax.legend(loc = 'upper right', fontsize=25)
# ax.set_xlim(0, 20)
# ax.set_ylim(60, 68)
plt.tight_layout()
plt.show()

################

Zc_vals = np.sqrt(Ll/(Cl + Cm_vals))

homogeneous_validity = (Zc_vals**2 - Zm_vals**2) / (Zc_vals**2 + Zm_vals**2)

plt.plot(d_vals * 1e6, homogeneous_validity, color = my_cmap3(1), marker = 'o', markersize = 8, linestyle = 'None')

ax = plt.gca()

# Customize the border settings
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2.5)
ax.spines['left'].set_linewidth(2.5)

# Adjust tick thickness and label size
ax.tick_params(axis='both', which='both', width=2.5)
ax.tick_params(axis='both', labelsize=25)

ax.set_xticks([0, 5, 10, 15, 20, 25, 30, 35, 40 ])
ax.set_xticklabels(['0', '5', '10', '15', '20', '25', '30', '35', '40'])

# Set axis labels and font size
ax.set_xlabel(r'$d$ (um)', fontsize=30)
ax.set_ylabel(r'$k_-/k_+$', fontsize=30)

# Add legend with fontsize 20

#ax.legend(loc = 'upper right', fontsize=25)
# ax.set_xlim(0, 20)
# ax.set_ylim(60, 68)
plt.tight_layout()
plt.show()