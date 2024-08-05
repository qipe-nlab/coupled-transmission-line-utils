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
from analytic_coupled_lines_calculator import calc_self_and_coupling_capacitance

import cap_util as cap

### param set 1: for high freq range

phase_vel = 1.2e8 #constants.c / np.sqrt(6.351)
Z0 = 65

d_val = 5e-6 # 12.5*1e-6 # separation between the coupled sections

Cm = cap.get_Cm(d_val)
Lm = cap.get_Lm(d_val)

d_vals = np.linspace(1, 25, 10) * 1e-6
Cm_vals = cap.get_Cm(d_vals) ### F per meter

# np.save('d_vals.npy', d_vals)
# np.save('Cm_vals.npy', Cm_vals)

C_val = cap.get_Cc(49e-6)
C_vals = cap.get_Cc(d_vals)

# plt.plot(d_vals * 1e6, (C_vals-Cm_vals) * 1e12, color = my_cmap2(7), linewidth = 7)
# plt.show()

L_val = cap.get_L(49e-6)
L_vals = cap.get_L(d_vals)

print('Cm:', Cm)
print('C_val:', C_val)

C_val2 = cap.get_Cc(d_val)
print('C_val2:', C_val2)

Z_val = np.sqrt(L_val/C_val)
Z_vals = np.sqrt(L_vals/C_vals)

# plt.plot(d_vals * 1e6, Z_vals, color = my_cmap2(7), linewidth = 7)
# plt.show()

v_val = 1/(np.sqrt(L_val*C_val))
v_vals = 1/np.sqrt(L_vals*C_vals)

print('Z_val:', Z_val)
print('v_val:', v_val)

#sys.exit()

w = 5e-6
s = 7.5e-6
eps_r = 11.7

cs, cm_vals_analytic = calc_self_and_coupling_capacitance(d_vals, w, s, eps_r)

# cs, cm_vals_analytic = calc_self_and_coupling_capacitance(d_vals, w, s, eps_r)

# d_vals_device = np.array([5.5, 5.5, 4.2, 3.8])*1e-6
# cs_, cm_for_device = calc_self_and_coupling_capacitance(d_vals_device, w, s, eps_r)
# print('cm_for_device:', cm_for_device)
#sys.exit()

#print('cm_vals_analytic:', cm_vals_analytic)

plt.figure(figsize=(8, 4))
plt.plot(d_vals * 1e6, Cm_vals * 1e12, color = my_cmap2(7), linewidth = 7, label = 'COMSOL')
plt.plot(d_vals * 1e6, cm_vals_analytic * 1e12, color = 'k', linewidth = 7, linestyle = '--', label = 'numeric')

### pF per meter - fF / mm

plt.xlabel(r'$d$ (um)', size = 45)
plt.ylabel(r'$c_m$ (fF/mm)', size = 45)
plt.legend(fontsize = 16)
#plt.title('Radiative T1 limit through detector line')

ax = plt.gca()

ax.tick_params(axis='both', which='major', labelsize=16)

# ax.set_yticks([1e-6,1e-3,1,1e3]) 
# ax.set_yticklabels(['1 ns','1 us','1 ms', '1 s'])
# ax.set_ylim([0.5e-6, 5e3])

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3)

plt.xticks([0, 5, 10, 15, 20, 25],[0, 5, 10, 15, 20, 25], size = 30)
plt.yticks([0, 5, 10, 15],['0', 5, 10, '15'], size = 30)
plt.grid(visible = True)
plt.xlim([0,25])
plt.ylim([0,15])
ax.tick_params(width=4)
plt.tight_layout()
plt.show()

# plt.close()

#sys.exit()

#######

print('cs:', cs)

plt.figure(figsize=(8, 4.25))
plt.plot(d_vals * 1e6, cm_vals_analytic/(cs[-1])*100, color = my_cmap2(7), linewidth = 7, label = 'Numerical') # (0,(5,4))
plt.plot(d_vals * 1e6, Cm_vals/C_vals[-1]*100, color = 'k', linewidth = 7, linestyle = (5,(5,5)), label = 'COMSOL')
#plt.plot(d_vals * 1e6, cm_vals_analytic/(cs), color = 'k', linewidth = 7, linestyle = '--')

plt.xlabel(r'$d$ (um)', size = 45)
plt.ylabel(r'$c_m/c~(\%)$', size = 45)

ax = plt.gca()

ax.tick_params(axis='both', which='major', labelsize=16)

# ax.set_yticks([1e-6,1e-3,1,1e3]) 
# ax.set_yticklabels(['1 ns','1 us','1 ms', '1 s'])
# ax.set_ylim([0.5e-6, 5e3])

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3)

plt.xticks([0, 5, 10, 15, 20, 25],[0, 5, 10, 15, 20, 25], size = 30)
plt.yticks(np.array([0, 0.025, 0.05, 0.075, 0.1, 0.125])*100, [0, '', 5, '', 10, ''], size = 30)
plt.grid(visible = True)
plt.xlim([0,25])
plt.ylim([0,12.5])
#plt.yscale('log')
plt.legend(fontsize = 20)
ax.tick_params(length = 5, width=2.5,direction='in')
plt.tight_layout()
plt.savefig('cm_c_ratio_ForPublication.pdf', format= 'pdf')
plt.show()

plt.plot(d_vals * 1e6, cs, color = my_cmap2(7), linewidth = 7)

plt.show()
sys.exit()

#######

phase_vel = 1.2e8 #constants.c / np.sqrt(6.351)
Z0 = 65

Cc = 1/(Z0*phase_vel)
Lc = Z0/phase_vel

d_val = 5e-6 # 12.5*1e-6 # separation between the coupled sections

Cm = cap.get_Cm(d_val)
Lm = cap.get_Lm(d_val)

print('c_m:', Cm)

d_vals = np.linspace(1, 25, 50) * 1e-6
Cm_vals = cap.get_Cm(d_vals) ### F per meter
Lm_vals = cap.get_Lm(d_vals) ### F per meter

L_ratios = Lm_vals/Lc
C_ratios = Cm_vals/Cc

Z_ms = (Lm_vals/Cm_vals)**0.5
Z_cs = (Lc/Cc)**0.5

Z_ratios_squared = (Z_ms/Z_cs)**2

C_vals =  cap.get_Cc(d_vals)
L_vals =  cap.get_L(d_vals)

Cm_vals =  cap.get_Cm(d_vals)
Lm_vals =  cap.get_Lm(d_vals)

######  split version

plt.figure(figsize=(9, 4))

plt.plot(d_vals * 1e6, L_ratios, color = my_cmap2(0), linewidth = 6, label = r'$l_m/l_c$')
plt.plot(d_vals * 1e6, C_ratios, color = my_cmap3(1), linewidth = 6, linestyle = ':', label = r'$c_m/c_c$')

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

plt.xticks([0, 5, 10, 15, 20, 25],[0, 5, 10, 15, 20, 25], size = 35)
#plt.yticks([0, 0.2, 0.4, 0.6, 0.8,1],[0, 0.2, 0.4, 0.6, 0.8,1], size = 35)
plt.yticks([0, 0.05, 0.1, 0.15],['0', '', '', '0.15'], size = 35)
plt.legend(loc='upper right', prop={'size': 35}, frameon = False)
ax.tick_params(width=3)
plt.tight_layout()
plt.xlim([0,25])
plt.ylim([0,0.15])
plt.grid(visible=True, axis='y')
plt.show()

###### Split two

plt.figure(figsize=(9, 4))

plt.plot(d_vals * 1e6, Z_ratios_squared, color = my_cmap3(3), linewidth = 6, label = r'$Z_m^2/Z_c^2$')
plt.plot(d_vals * 1e6, Z_vals/Z_val, color = my_cmap2(0), linewidth = 6, label = r'$Z_c/Z_0$')
plt.plot(d_vals * 1e6, v_vals/v_val, color = my_cmap3(1), linewidth = 6, linestyle = ':', label = r'$v_c/v$')

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

plt.xticks([0, 5, 10, 15, 20, 25],[0, 5, 10, 15, 20, 25], size = 35)
#plt.yticks([0, 0.2, 0.4, 0.6, 0.8,1],[0, 0.2, 0.4, 0.6, 0.8,1], size = 35)
plt.yticks([0.95, 1, 1.05, 1.1],['0.95', '1', '1.05', '1.1'], size = 35)
plt.legend(loc='upper right', prop={'size': 35}, frameon = False)
ax.tick_params(width=3)
plt.tight_layout()
plt.xlim([0,25])
plt.ylim([0.95,1.1])
plt.grid(visible=True, axis='y')
plt.show()

#######