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
import cap_util as cap

### param set 1: for high freq range

phase_vel = 1.2e8 #constants.c / np.sqrt(6.351)
Z0 = 65

d_val = 5e-6 # 12.5*1e-6 # separation between the coupled sections

Cm = cap.get_Cm(d_val)
Lm = cap.get_Lm(d_val)

d_vals = np.linspace(1, 25, 50) * 1e-6
Cm_vals = cap.get_Cm(d_vals) ### F per meter

C_val = cap.get_Cc_plus_Cm(49e-6)
C_vals = cap.get_Cc_plus_Cm(d_vals)

# plt.plot(d_vals * 1e6, (C_vals-Cm_vals) * 1e12, color = my_cmap2(7), linewidth = 7)
# plt.show()

L_val = cap.get_L(49e-6)
L_vals = cap.get_L(d_vals)

print('Cm:', Cm)
print('C_val:', C_val)

C_val2 = cap.get_Cc_plus_Cm(d_val)
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

plt.figure(figsize=(8, 4))
plt.plot(d_vals * 1e6, Cm_vals * 1e12, color = my_cmap2(7), linewidth = 7)

### pF per meter - fF / mm

plt.xlabel(r'$d$ (um)', size = 45)
plt.ylabel(r'$c_m$ (fF/mm)', size = 45)
#plt.legend(fontsize = 16)
#plt.title('Radiative T1 limit through detector line')

ax = plt.gca()

ax.tick_params(axis='both', which='major', labelsize=16)

# ax.set_yticks([1e-6,1e-3,1,1e3]) 
# ax.set_yticklabels(['1 ns','1 us','1 ms', '1 s'])
# ax.set_ylim([0.5e-6, 5e3])

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3)

plt.xticks([0, 5, 10, 15, 20, 25],[0, 5, 10, 15, 20, 25], size = 45)
plt.yticks([0, 5, 10, 15],['0', '', '', '15'], size = 45)
plt.grid(visible = True)
plt.xlim([0,25])
plt.ylim([0,15])
ax.tick_params(width=4)
plt.tight_layout()
plt.show()

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

plt.figure(figsize=(6, 5))

L_ratios = Lm_vals/Lc
C_ratios = Cm_vals/Cc

Z_ms = (Lm_vals/Cm_vals)**0.5
Z_cs = (Lc/Cc)**0.5

Z_ratios = Z_ms/Z_cs

plt.plot(d_vals * 1e6, L_ratios, color = my_cmap2(0), linewidth = 6, label = r'$l_m/l_c$')
plt.plot(d_vals * 1e6, C_ratios, color = my_cmap3(1), linewidth = 6, linestyle = ':', label = r'$c_m/c_c$')
plt.plot(d_vals * 1e6, Z_ratios, color = my_cmap2(6), linewidth = 6, label = r'$Z_m/Z_c$')
plt.plot(d_vals * 1e6, Z_vals/Z_val, color = my_cmap3(3), linewidth = 6, label = r'$Z_c/Z_0$')
plt.plot(d_vals * 1e6, v_vals/v_val, color = my_cmap3(0), linewidth = 6, linestyle = ':', label = r'$v_c/v$')

### pF per meter - fF / mm

C_vals =  cap.get_C(d_vals)
L_vals =  cap.get_L(d_vals)

Cm_vals
Lm_vals =  cap.get_Lm(d_vals)

test = L_vals*C_vals - Cm_vals*Lm_vals

test = test*((3*10**8)/(3/1.2))**2

plt.plot(d_vals, test)
plt.show()

###

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

plt.xticks([1, 5, 10, 15, 20, 25],[1, 5, 10, 15, 20, 25], size = 35)
#plt.yticks([0, 0.2, 0.4, 0.6, 0.8,1],[0, 0.2, 0.4, 0.6, 0.8,1], size = 35)
plt.yticks([0, 1],[0, 1], size = 35)
plt.legend(loc='center right', prop={'size': 35}, frameon = False)
ax.tick_params(width=3)
plt.tight_layout()
plt.xlim([1,25])
plt.ylim([0,1.1])
plt.grid(visible=True, axis='y')
plt.show()

###### test split version

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

plt.plot(d_vals * 1e6, Z_ratios, color = my_cmap3(3), linewidth = 6, label = r'$Z_m/Z_c$')
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

print('Cm:', Cm)
print('Lm:', Lm)

# l_Rf = 1.7e-3
# l_Rn = 0.96e-3
# l_Gf = 1.2e-3
# l_Gn = 1.5e-3
# l_c = 0.25e-3

l_Rf = 1.7e-3
l_Rn = 0.96e-3
l_Gf = 1.7e-3
l_Gn = 1.0e-3
l_c = 0.275e-3

L_readout = l_Gn + l_c + l_Gf

readout_omega = lambda_quarter_omega(L_readout, phase_vel=phase_vel)

L_purcell = l_Rn + l_c + l_Rf

purcell_omega = lambda_quarter_omega(L_purcell, phase_vel=phase_vel)

print('readout_omega:', readout_omega / (2*np.pi*1e9))
print('purcell_omega:', purcell_omega / (2*np.pi*1e9))

print('Lm:', Lm)
print('Cm:', Cm)
print('l_Rf:', l_Rf)
print('l_Rn:', l_Rn)
print('l_Gf:', l_Gf)
print('l_Gn:', l_Gn)
print('l_c:', l_c)
#sys.exit()

test_notch_freq = find_notch_filter_frequency(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=phase_vel, Z0=Z0)

test_notch_freq_analytic = find_notch_filter_frequency_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=phase_vel, Z0=Z0)

test_notch_freq_rule_of_thumb = notch_filter_frequency_rule_of_thumb(l_c, l_Gf,l_Rf, Cm, phase_vel=phase_vel, Z0=Z0)

print('rule_of_thumb_freq1 (GHz):', test_notch_freq_rule_of_thumb/ (2*np.pi*1e9))

test_omega_1 = lambda_quarter_omega(l_c + l_Gf + l_Gn, phase_vel=phase_vel)

test_omega_2 = lambda_quarter_omega(l_c + l_Rf + l_Rn, phase_vel=phase_vel)

J_val = J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm)

# print('test_notch_freq (GHz):', test_notch_freq/(2*np.pi*1e9))
print('test_omega_1 (GHz):', test_omega_1/(2*np.pi*1e9))
print('test_omega_2 (GHz):', test_omega_2/(2*np.pi*1e9))
print('J_val (MHz):', J_val/(2*np.pi*1e6))

omegas = np.arange(7.5115, 11.0115, 0.005) * 1e9 * 2*np.pi

test_Z_transfer_exact = Z_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

# # test_S_transfer_sym_3_lines_exact = S_transfer_sym_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

# # print('test_S_transfer_sym_3_lines_exact:', test_S_transfer_sym_3_lines_exact)

# # S_transfer = test_S_transfer_sym_3_lines_exact[:, 0, 1]

# # print('S_transfer:', S_transfer)

# # plt.plot(omegas/(2*np.pi*1e9), np.abs(S_transfer))
# # plt.yscale('log')
# # plt.show()

# Z_test_exact_vals = voltage_at_source_location_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, Z0, phase_vel, omegas)
# Z_test_approx_vals = Z_trans_along_shorted_tl(Z0, phase_vel, l_Gn+l_Gf+l_c, l_Gn, omegas)

# plt.plot(omegas/(2*np.pi * 1e9), np.imag(Z_test_exact_vals), color = my_cmap3(1), label = 'Exact solution', linewidth = 3)
# plt.plot(omegas/(2*np.pi * 1e9), np.imag(Z_test_approx_vals), color = my_cmap2(7), label = 'Single transmission line solution', linewidth = 3, alpha = 0.6)

# plt.yscale('log')

# plt.xlabel('Frequency (GHz)', size = 20)
# plt.ylabel(r'$V_1 / I_{in}$' + r' $(\Omega)$', size = 20)
# plt.legend(fontsize = 16)
# #plt.title('Radiative T1 limit through detector line')

# ax = plt.gca()

# ax.tick_params(axis='both', which='major', labelsize=16)

# # ax.set_yticks([1e-6,1e-3,1,1e3]) 
# # ax.set_yticklabels(['1 ns','1 us','1 ms', '1 s'])
# # ax.set_ylim([0.5e-6, 5e3])

# for axis in ['top','bottom','left','right']:
#     ax.spines[axis].set_linewidth(2)

# ax.tick_params(width=3)
# plt.tight_layout()
# plt.show()

###

# Z11_exact = Z_input_3_lines_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

# Z11_LE_circuit = Z_input_equivalent_LE_circuit(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

# Z11_LE_circuit_symbolic = Z_input_equivalent_LE_circuit(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

# plt.plot(omegas/(2*np.pi * 1e9), np.abs(Z11_exact), color = my_cmap3(1), label = 'Distributed circuit', linewidth = 3)
# #plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_weak_coupling), color = 'r', linestyle = '--', label = 'Z transfer weak coupling model')
# plt.plot(omegas/(2*np.pi*1e9), np.abs(Z11_LE_circuit), color = my_cmap2(7), label = 'Equivalent circuit', linewidth = 3, alpha = 0.85)
# #plt.plot(omegas/(2*np.pi*1e9), np.abs(Z11_LE_circuit_symbolic), color = my_cmap3(7), linestyle = '--', label = 'Equivalent circuit - symbolic', linewidth = 3, alpha = 0.5)

plt.yscale('log')
plt.legend(loc = 'upper left', fontsize = 18, frameon=False)
plt.xlabel('Frequency (GHz)', size = 20)
plt.ylabel(r'Impedance $Z_{11}$ ($\Omega$)', size = 20)
#plt.title('Z transfer function for different models')
# Customize the border settings
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(3)
ax.spines['left'].set_linewidth(3)

# ax.set_xticks([8, 9, 10, 11])

# #ax.set_yticks([1e-2, 1, 100, 1e4])

# # Adjust tick thickness and label size
# ax.tick_params(axis='both', which='both', width=2.5)
# ax.tick_params(axis='both', labelsize=20)
# plt.tight_layout()

# plt.show()

###

## Z21!

plt.figure(figsize=(7.5, 4.25))

Z_test_exact_vals = voltage_at_source_location_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, Z0, phase_vel, test_notch_freq_rule_of_thumb)
Z_test_approx_vals = Z_trans_along_shorted_tl(Z0, phase_vel, l_Gn+l_Gf+l_c, l_Gn, test_notch_freq_rule_of_thumb)

Z_ratios = Z_test_exact_vals/Z_test_approx_vals

#test_Z_transfer_weak_coupling = Z_transfer_weak_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

#test_Z_equiv_LE_circuit = Z_transfer_equivalent_LE_circuit(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

test_Z_equiv_LE_circuit_symbolic = Z_transfer_equivalent_LE_circuit_from_symbolic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Cm, omegas, phase_vel=phase_vel, Z0=Z0)

# my_cmap2(7)

plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_exact), color = my_cmap3(1), label = 'Distributed circuit', linewidth = 3, zorder = 2)
#plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_weak_coupling), color = 'r', linestyle = '--', label = 'Z transfer weak coupling model')
#plt.plot(omegas/(2*np.pi*1e9), np.abs(test_Z_equiv_LE_circuit), color = my_cmap2(7), label = 'Equivalent circuit', linewidth = 3, alpha = 0.85)
plt.plot(omegas/(2*np.pi*1e9), np.abs(test_Z_equiv_LE_circuit_symbolic), color = my_cmap2(7), label = 'Equivalent circuit', linestyle = '--', linewidth = 3, alpha = 0.7, zorder = 11)

#print('test_notch_freq_rule_of_thumb/(2*np.pi*1e9):', test_notch_freq_rule_of_thumb/(2*np.pi*1e9))
#plt.vlines(test_notch_freq_analytic/(2*np.pi*1e9), 0.008, 2, color = my_cmap(7), linestyle = 'dotted', linewidth = 3)
#plt.vlines(test_notch_freq_rule_of_thumb/(2*np.pi*1e9), 0.008, 2, color = my_cmap(4), linestyle = 'dotted', linewidth = 3, label = 'Notch eq.',zorder=10)

### add direct cap result

print('L_readout:', L_readout)
print('L_purcell:', L_purcell)

l_Rf = 1.725e-3
l_Rn = 1.21e-3
l_Gf = 1.225e-3
l_Gn = 1.75e-3

print('l_Rf + l_Rn:', l_Rf + l_Rn)
print('l_Gf + l_Gn:', l_Gf + l_Gn)

C_direct_val = get_eff_Cf_from_J(0, l_Gf, l_Gn, l_Rf, l_Rn, J_val, phase_vel=3*10**8/2.5, Z0=65)

C_direct_val = 2.15*C_direct_val

J_cap_val = J_coupling_direct_cap(l_Gf, l_Gn, l_Rf, l_Rn, C_direct_val, phase_vel=3*10**8/2.5, Z0=65)

print('C_direct_val:', C_direct_val*1e15)
print('J_cap_val (MHz):', J_cap_val/(2*np.pi*1e6))

test_Z_transfer_exact = Z_transfer_direct_cap_exact(l_Gf, l_Gn, l_Rf, l_Rn, C_direct_val, omegas, phase_vel=phase_vel, Z0=Z0)

test_Z_equiv_LE_circuit_symbolic = Z_transfer_direct_cap_equivalent_LE_circuit(l_Gf, l_Gn, l_Rf, l_Rn, C_direct_val, omegas, phase_vel=phase_vel, Z0=Z0)

# my_cmap2(7)

plt.plot(omegas/(2*np.pi * 1e9), np.abs(test_Z_transfer_exact), color = my_cmap3(1), label = 'cap', linewidth = 4, zorder = 1)
plt.plot(omegas/(2*np.pi*1e9), np.abs(test_Z_equiv_LE_circuit_symbolic), color = 'lightskyblue', label = 'Equiv. cap', linestyle = '--',  linewidth = 4, alpha = 0.7, zorder = 10)

plt.yscale('log')
plt.legend(loc = 'upper left', fontsize = 22.5, frameon = False)
plt.xlabel('Frequency (GHz)', size = 25)
plt.ylabel(r'Impedance $Z_{21}$ ($\Omega$)', size = 25)
plt.legend(loc = 'upper left', fontsize = 18, frameon=False)
plt.xlabel('Frequency (GHz)', size = 20)
plt.ylabel(r'Impedance $Z_{21}$ ($\Omega$)', size = 20)
#plt.title('Z transfer function for different models')
# Customize the border settings
ax = plt.gca()
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
ax.spines['top'].set_linewidth(3)
ax.spines['right'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.spines['left'].set_linewidth(3)

ax.set_xticks([8, 9, 10, 11])

ax.set_yticks([1e-2, 1, 100, 1e4])
#ax2.set_yticklabels(['0', '5', '10', '15'])

# Adjust tick thickness and label size
ax.tick_params(axis='both', which='both', width=2.5)
ax.tick_params(axis='both', labelsize=25)
#plt.grid(True)
plt.tight_layout()

plt.show()

#sys.exit()

C_q = 50e-15

### param set 1: for high freq range
C_g = 5.5e-15 #3.65e-15
C_ext = 10e-15 #24e-15

### param set 2: for low freq range
# C_g = 6e-15
# C_ext = 52.5e-15

omega_q = test_notch_freq 

test_g_val = g_coupling(omega_q, C_q, C_g, l_c, l_Gf, l_Gn, phase_vel=3*10**8/2.5, Z0=65)
test_J_val = J_coupling(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=3*10**8/2.5, Z0=65)
test_J_val_symbolic = J_coupling_analytic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, Cm, phase_vel=3*10**8/2.5, Z0=65)

test_k_val = k_readout(C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)

print('test_g_val (MHz):', test_g_val/(2*np.pi*1e6))
print('test_J_val (MHz):', test_J_val/(2*np.pi*1e6))
print('test_J_val_symbolic (MHz):', test_J_val_symbolic/(2*np.pi*1e6))

J_ratios = test_J_val/test_J_val_symbolic

print('J_ratios:', J_ratios)
print('Z_ratios:', Z_ratios)

print('test_k_val (MHz):', test_k_val/(2*np.pi*1e6))

omegas = np.arange(7, 11, 0.005) * 1e9 * 2*np.pi

test_T1_radiative_exact = qubit_radiative_decay_sym_3_lines_exact(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)
test_T1_radiative_equivalent_LE_circuit = qubit_radiative_decay_equivalent_LE_circuit(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)
test_T1_radiative_equivalent_LE_circuit_without_notch = qubit_radiative_decay_equivalent_LE_circuit_without_notch(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)
test_T1_radiative_equivalent_LE_circuit_single_resonator = qubit_radiative_decay_equivalent_LE_circuit_single_resonator(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)

test_T1_radiative_equivalent_LE_circuit_approx = qubit_radiative_decay_equivalent_LE_circuit_approximate(C_q, C_g, C_ext, l_c, l_Gf, l_Gn, l_Rf, l_Rn, Lm, Cm, omegas, phase_vel=3*10**8/2.5, Z0=65, Zline = 50)

plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_exact * 1e3, color = my_cmap3(1), linewidth = 3, label = 'with intrinsic notch - exact')
plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit * 1e3, color = my_cmap2(7), linewidth = 3, label = 'with intrinsic notch - e.c.', alpha = 0.5)
plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit_approx * 1e3, color = my_cmap2(7), linewidth = 3, linestyle = '-.', label = 'with intrinsic notch - e.c. approx', alpha = 0.5)
plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit_without_notch * 1e3, color = my_cmap2(7), linewidth = 3, linestyle = '--', label = 'without intrinsic notch - e.c.', alpha = 0.5)
# plt.plot(omegas/(2*np.pi * 1e9), test_T1_radiative_equivalent_LE_circuit_single_resonator * 1e3, color = 'b', linestyle = '--', label = 'equiv. single resonator circuit')

plt.yscale('log')
plt.xlabel('Frequency (GHz)', size = 20)
plt.ylabel(r'Purcell decay $T_1$ limit', size = 20)
plt.legend(fontsize = 16)
#plt.title('Radiative T1 limit through detector line')

ax = plt.gca()

ax.tick_params(axis='both', which='major', labelsize=16)

ax.set_yticks([1e-6,1e-3,1,1e3]) 
ax.set_yticklabels(['1 ns','1 us','1 ms', '1 s'])
ax.set_ylim([0.5e-6, 5e3])

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)

ax.tick_params(width=3)
plt.tight_layout()

plt.show()

enhancement_factor = test_T1_radiative_exact/test_T1_radiative_equivalent_LE_circuit_without_notch
enhancement_factor_approx = test_T1_radiative_equivalent_LE_circuit/test_T1_radiative_equivalent_LE_circuit_without_notch

test_enhancement_bandwidth = notch_enhancement_bandwidth(l_c, l_Gf, l_Gn, l_Rf, l_Rn, T1_enhancement_fact = 10, phase_vel = 3*10**8/2.5)

print('test_enhancement_bandwidth (GHz):', test_enhancement_bandwidth / (2*np.pi*1e9))

test_enhancement_bandwidth_exact = notch_enhancement_bandwidth_exact(l_c, l_Gf, l_Gn, l_Rf, l_Rn, T1_enhancement_fact = 10, phase_vel = 3*10**8/2.5)

print('test_enhancement_bandwidth_exact (GHz):', test_enhancement_bandwidth_exact / (2*np.pi*1e9))

plt.plot(omegas/(2*np.pi * 1e9), enhancement_factor, color = my_cmap2(7), linewidth = 3, label = 'exact solution')
#plt.plot(omegas/(2*np.pi * 1e9), enhancement_factor_approx, color = my_cmap2(7), linewidth = 3, label = 'equivalent circuit solution', alpha = 0.5)

omega_r = lambda_quarter_omega(l_c + l_Gf + l_Gn, phase_vel=phase_vel)
omega_p = lambda_quarter_omega(l_c + l_Rf + l_Rn, phase_vel=phase_vel)

predicted_notch_enhancement = enhancement_factor_symbolic(l_c, l_Gf, l_Gn, l_Rf, l_Rn, omegas)
#plt.plot(omegas/(2*np.pi * 1e9), predicted_notch_enhancement, color = 'g')

exact_bandwidth = find_notch_enhancement_bandwidth(omegas, enhancement_factor, 10)

print('exact_bandwidth (GHz):', exact_bandwidth / (2*np.pi*1e9))

plt.yscale('log')
plt.xlabel('Frequency (GHz)', size = 20)
plt.ylabel(r'Purcell $T_1$ enhancement factor', size = 20)
#plt.legend(fontsize = 16)

ax = plt.gca()

ax.tick_params(axis='both', which='major', labelsize=16)

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)

ax.tick_params(width=3)
plt.tight_layout()
plt.xlim(7, 9.75)
plt.hlines(10, 6, 11, color = 'k', linestyle = '--')
plt.yscale('log')
plt.yticks([1, 10, 100, 1000, 10**4, 10**5], ['1', '10', '100', r'$10^3$', r'$10^4$', r'$10^5$'])
plt.ylim(0.5, 5*10**5)
plt.show()

##### test predicted T1 limits from circuit model:
