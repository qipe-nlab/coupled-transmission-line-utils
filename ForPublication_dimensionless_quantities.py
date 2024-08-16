import seaborn as sns
from scipy import constants

hex_codes = sns.color_palette().as_hex()
hex_codes2 = sns.color_palette('husl', 9).as_hex()

from matplotlib.colors import ListedColormap

#flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
my_cmap = ListedColormap(hex_codes)
my_cmap2 = ListedColormap(hex_codes2)
my_cmap3 = ListedColormap(["#9b59b6", "#34495e", "#3498db", "#95a5a6", "#e74c3c"])

### example solutions

from exact_coupled_transmission_line_eqn_solver import *
from analytic_coupled_lines_calculator import calc_self_and_coupling_capacitance, calc_self_and_coupling_inductance

import cap_util as cap

### param set 1: for high freq range

w = 5e-6
s = 7.5e-6
eps_r = 11.7

show_conformal_mapping_result = False

### conformal mapping result

if show_conformal_mapping_result:

    d_vals = np.linspace(1, 25, 10) * 1e-6

    large_dist = 500*1e-6
    cs_large_dist, cm_large_dist = calc_self_and_coupling_capacitance(large_dist, w, s, eps_r)
    ls_large_dist, lm_large_dist = calc_self_and_coupling_inductance(large_dist, w, s)

    Z0 = np.sqrt(ls_large_dist/cs_large_dist)
    phase_vel = np.sqrt(1/(ls_large_dist * cs_large_dist))

    Cc_vals, Cm_vals = calc_self_and_coupling_capacitance(d_vals, w, s, eps_r) ### F per meter
    Lc_vals, Lm_vals =  calc_self_and_coupling_inductance(d_vals, w, s)

    # print('Cc_vals:', Cc_vals)
    # print('Cm_vals:', Cm_vals)
    # print('Lc_vals:', Lc_vals)
    # print('Lm_vals:', Lm_vals)

    Cc_vals = np.array(Cc_vals)
    Cm_vals = np.array(Cm_vals)
    Lc_vals = np.array(Lc_vals)
    Lm_vals = np.array(Lm_vals)

    Zc_vals = np.sqrt(Lc_vals/(Cc_vals+Cm_vals))
    vc_vals = 1/np.sqrt(Lc_vals*(Cc_vals+Cm_vals))

    L_ratios = Lm_vals/Lc_vals
    C_ratios = Cm_vals/Cc_vals

    Z_ms = (Lm_vals/Cm_vals)**0.5
    Z_cs = (Lc_vals/(Cc_vals+Cm_vals))**0.5

    k = (1-(Lm_vals/Lc_vals)**2)**0.5

    plt.figure(figsize=(10, 4))

    plt.plot(d_vals * 1e6, (Z_ms/Z_cs)**2, color = my_cmap3(3), linewidth = 6, linestyle = ':', label = r'$Z_m^2/Z_c^2$')
    plt.plot(d_vals * 1e6, Zc_vals/Z0, color = my_cmap2(0), linewidth = 6, linestyle = ':', label = r'$Z_c/Z_0$')
    plt.plot(d_vals * 1e6, vc_vals/phase_vel, color = my_cmap3(1), linewidth = 6, linestyle = ':', label = r'$v_c/v$')
    plt.plot(d_vals * 1e6, k, color = my_cmap3(5), linewidth = 6, linestyle = ':', label = r'$k$')

    plt.xlabel(r'$d$ (um)', size = 35)
    #plt.ylabel(r'$c_m$ (fF/mm)', size = 35)
    #plt.legend(fontsize = 16)
    #plt.title('Radiative T1 limit through detector line')

    ax = plt.gca()

    ax.tick_params(axis='both', which='major', labelsize=16)

    # ax.set_yticks([1e-6,1e-3,1,1e3]) 
    # ax.set_yticklabels(['1 ns','1 us','1 ms', '1 s'])
    # ax.set_ylim([0.5e-6, 5e3])

    # for axis in ['top','bottom','left','right']:
    #     ax.spines[axis].set_linewidth(3)

    # plt.xticks([0, 5, 10, 15, 20, 25],[0, 5, 10, 15, 20, 25], size = 35)
    # #plt.yticks([0, 0.2, 0.4, 0.6, 0.8,1],[0, 0.2, 0.4, 0.6, 0.8,1], size = 35)
    # plt.yticks([0.95, 1, 1.05, 1.1],['0.95', '1', '1.05', '1.1'], size = 35)
    # plt.legend(loc='upper right', prop={'size': 35}, frameon = False)
    # ax.tick_params(width=3)
    # plt.tight_layout()
    # plt.xlim([0,25])
    # plt.ylim([0.95,1.05])
    # plt.grid(visible=True, axis='y')
    # plt.axvline(x=3.8, color='k', linestyle='-')
    # plt.axvline(x=5.5, color='k', linestyle='-')
    #plt.savefig('justifying_CPW_assumptions_test.pdf', format = 'pdf')
    #plt.show()

### COMSOL result

d_vals = np.linspace(1, 25, 50) * 1e-6
Cm_vals = cap.get_Cm(d_vals) ### F per meter
Lm_vals =  cap.get_Lm(d_vals)

Cc_vals = cap.get_Cc(d_vals)
Lc_vals = cap.get_L(d_vals)

Z0 = np.sqrt(Lc_vals[-1]/(Cc_vals[-1]+Cm_vals[-1]))
phase_vel = np.sqrt(1/(Lc_vals[-1] * (Cc_vals[-1]+Cm_vals[-1])))

# #######

Zc_vals = np.sqrt(Lc_vals/(Cc_vals+Cm_vals))
vc_vals = 1/np.sqrt(Lc_vals*(Cc_vals+Cm_vals))

L_ratios = Lm_vals/Lc_vals
C_ratios = Cm_vals/Cc_vals

Z_ms = (Lm_vals/Cm_vals)**0.5
Z_cs = (Lc_vals/(Cc_vals+Cm_vals))**0.5

k = (1-(Lm_vals/Lc_vals)**2)**0.5

### plot

if not show_conformal_mapping_result:

    plt.figure(figsize=(10, 3.5))

#plt.figure(figsize=(10, 4))

# plt.plot(d_vals * 1e6, Cm_vals/Cc_vals, color = my_cmap3(3), linewidth = 5, linestyle = '-', label = r'$Z_m^2/Z_c^2$', alpha = 0.75)
# plt.plot(d_vals * 1e6, Lm_vals/Lc_vals, color = my_cmap3(3), linewidth = 5, linestyle = '-', label = r'$Z_m^2/Z_c^2$', alpha = 0.75)

plt.plot(d_vals * 1e6, (Z_ms/Z_cs), color = my_cmap3(3), linewidth = 5, linestyle = '-', label = r'$Z_m/Z_c$', alpha = 0.75)
plt.plot(d_vals * 1e6, Zc_vals/Z0, color = my_cmap2(0), linewidth = 5, linestyle = '-', label = r'$Z_c/Z_0$', alpha = 0.75)
plt.plot(d_vals * 1e6, vc_vals/phase_vel, color = my_cmap3(1), linewidth = 5, linestyle = '-', label = r'$v_c/v$', alpha = 0.75)
plt.plot(d_vals * 1e6, k, color = my_cmap3(5), linewidth = 5, linestyle = '-', label = r'$k$', alpha = 0.75)

plt.xlabel(r'$d$ (um)', size = 35)
#plt.ylabel(r'$c_m$ (fF/mm)', size = 35)
#plt.legend(fontsize = 16)
#plt.title('Radiative T1 limit through detector line')

ax = plt.gca()

ax.tick_params(axis='both', which='major', labelsize=16)

# ax.set_yticks([1e-6,1e-3,1,1e3]) 
# ax.set_yticklabels(['1 ns','1 us','1 ms', '1 s'])
# ax.set_ylim([0.5e-6, 5e3])

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3)

plt.xticks([0, 5, 10, 15, 20, 25],[0, 5, 10, 15, 20, 25], size = 25)
#plt.yticks([0, 0.2, 0.4, 0.6, 0.8,1],[0, 0.2, 0.4, 0.6, 0.8,1], size = 35)
plt.yticks([0.98, 0.99, 1, 1.01, 1.02],['0.98', '0.99', '1', '1.01', '1.02'], size = 25)
plt.legend(loc='upper right', prop={'size': 25}, frameon = False)
ax.tick_params(width=3, length = 8, direction = 'in')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
plt.tight_layout()
plt.xlim([0,25])
plt.ylim([0.99,1.02])
plt.grid(visible=True, axis='y')
plt.axvline(x=3.8, color='k', linestyle=':')
plt.axvline(x=5.5, color='k', linestyle=':')
plt.savefig('justifying_CPW_assumptions_test.pdf', format = 'pdf')
plt.show()

#######