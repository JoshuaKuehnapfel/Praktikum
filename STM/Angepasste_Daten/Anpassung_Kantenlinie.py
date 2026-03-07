
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams['font.size'] = 24.0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelsize'] = 'medium'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['lines.linewidth'] = 2.0


H_200 = np.loadtxt('C:/VS_Studio_Latex/Uni (lLokale Kopie)/Praktikum/STM_angepasst/HOPG_200nm_Profillinie.txt', skiprows=4, delimiter=';', usecols=(0, 1))
H_23 = np.loadtxt('C:/VS_Studio_Latex/Uni (lLokale Kopie)/Praktikum/STM_angepasst/HOPG_23nm_Profillinie.txt', skiprows=4, delimiter=';', usecols=(0, 1))
H_12 = np.loadtxt('C:/VS_Studio_Latex/Uni (lLokale Kopie)/Praktikum/STM_angepasst/HOPG_12nm_Profillinie.txt', skiprows=4, delimiter=';', usecols=(0, 1))

def line(x, m, b):
    return m * x + b

def eval_line_and_var(x0, popt, pcov):
    m, b = popt
    val = m * x0 + b
    var = (x0 ** 2) * pcov[0, 0] + 2 * x0 * pcov[0, 1] + pcov[1, 1]
    return val, var

edge = np.where(H_200[:, 0] > 1.0e-8)[0][0]
popt_200_l, pcov_200_l = curve_fit(line, H_200[:edge-2, 0], H_200[:edge-2, 1])
popt_200_r, pcov_200_r = curve_fit(line, H_200[edge:, 0], H_200[edge:, 1])
val_l_200, var_l_200 = eval_line_and_var(1.0e-8, popt_200_l, pcov_200_l)
val_r_200, var_r_200 = eval_line_and_var(1.0e-8, popt_200_r, pcov_200_r)
Hoehe200 = val_l_200 - val_r_200
Hoehe200_var = var_l_200 + var_r_200
Hoehe200_err = np.sqrt(Hoehe200_var)
mean_H_200 = np.mean(H_200[:, 1])
x_werte_l = np.linspace(min(H_200[:edge-2, 0]), 1.0e-8, 50)
x_werte_r = np.linspace(1.0e-8, max(H_200[edge:, 0]), 50)
plt.figure(figsize=(11, 6.5))
plt.plot(H_200[:, 0], H_200[:, 1]-mean_H_200, label='200nm', marker='x', linestyle='dashed')
plt.plot(x_werte_l, line(x_werte_l, *popt_200_l)-mean_H_200, label='200nm links')
plt.plot(x_werte_r, line(x_werte_r, *popt_200_r)-mean_H_200, label='200nm rechts')
plt.vlines(1.0e-8, min(H_200[:, 1]-mean_H_200), max(H_200[:, 1]-mean_H_200), colors='gray', linestyles='dashed', label='Kantenlinie')
plt.xlabel('x (m)')
plt.title(f'Kantenanpassung 200 nm — Höhe = {Hoehe200*1e9:.2f} ± {Hoehe200_err*1e9:.2f} nm')
plt.ylabel('z (m)')
plt.grid(True)
plt.legend()
plt.savefig('Kantenanpassung_200nm.svg', dpi=300)

edge = np.where(H_23[:,0] > 0.174e-8)[0][0]
popt_23_l, pcov_23_l = curve_fit(line, H_23[:edge, 0], H_23[:edge, 1])
popt_23_r, pcov_23_r = curve_fit(line, H_23[edge+1:, 0], H_23[edge+1:, 1])
val_l_23, var_l_23 = eval_line_and_var(0.174e-8, popt_23_l, pcov_23_l)
val_r_23, var_r_23 = eval_line_and_var(0.174e-8, popt_23_r, pcov_23_r)
Hoehe23 = val_l_23 - val_r_23
Hoehe23_var = var_l_23 + var_r_23
Hoehe23_err = np.sqrt(Hoehe23_var)
mean_H_23 = np.mean(H_23[:, 1])
x_werte_l = np.linspace(min(H_23[:edge, 0]), 0.174e-8, 50)
x_werte_r = np.linspace(0.174e-8, max(H_23[edge+1:, 0]), 50)

plt.figure(figsize=(11, 6.5))
plt.plot(H_23[:, 0], H_23[:, 1]-mean_H_23, label='23nm', marker='x', linestyle='dashed')
plt.plot(x_werte_l, line(x_werte_l, *popt_23_l)-mean_H_23, label='23nm links')
plt.plot(x_werte_r, line(x_werte_r, *popt_23_r)-mean_H_23, label='23nm rechts')
plt.vlines(0.174e-8, min(H_23[:, 1]-mean_H_23), max(H_23[:, 1]-mean_H_23), colors='gray', linestyles='dashed', label='Kantenlinie')
plt.xlabel('x (m)')
plt.ylabel('z (m)')
plt.title(f'Kantenanpassung 23 nm — Höhe = {Hoehe23*1e9:.2f} ± {Hoehe23_err*1e9:.2f} nm')
plt.legend()
plt.grid(True)
plt.savefig('Kantenanpassung_23nm.svg', dpi=300)

edge = np.where(H_12[:, 0] > 0.136e-8)[0][0]
popt_12_l, pcov_12_l = curve_fit(line, H_12[15:edge-1, 0], H_12[15:edge-1, 1])
popt_12_r, pcov_12_r = curve_fit(line, H_12[edge:, 0], H_12[edge:, 1])
val_l_12, var_l_12 = eval_line_and_var(0.136e-8, popt_12_l, pcov_12_l)
val_r_12, var_r_12 = eval_line_and_var(0.136e-8, popt_12_r, pcov_12_r)
Hoehe12 = val_l_12 - val_r_12
Hoehe12_var = var_l_12 + var_r_12
Hoehe12_err = np.sqrt(Hoehe12_var)
mean_H_12 = np.mean(H_12[:, 1])
x_werte_l = np.linspace(min(H_12[15:edge-1, 0]), 0.136e-8, 50)
x_werte_r = np.linspace(0.136e-8, max(H_12[edge:, 0]), 50)

plt.figure(figsize=(11, 6.5))
plt.plot(H_12[15:, 0], H_12[15:, 1]-mean_H_12, label='12nm', marker='x', linestyle='dashed')
plt.plot(x_werte_l, line(x_werte_l, *popt_12_l)-mean_H_12, label='12nm links')
plt.plot(x_werte_r, line(x_werte_r, *popt_12_r)-mean_H_12, label='12nm rechts')
plt.vlines(0.136e-8, min(H_12[:, 1]-mean_H_12), max(H_12[:, 1]-mean_H_12), colors='gray', linestyles='dashed', label='Kantenlinie')
plt.xlabel('x (m)')
plt.ylabel('z (m)')
plt.title(f'Kantenanpassung 12 nm — Höhe = {Hoehe12*1e9:.2f} ± {Hoehe12_err*1e9:.2f} nm')
plt.legend()
plt.grid(True)
plt.savefig('Kantenanpassung_12nm.svg', dpi=300)

print(f'Höhe 200nm: {Hoehe200*1e9:.2f} ± {Hoehe200_err*1e9:.2f} nm')
print(f'Höhe 23nm: {Hoehe23*1e9:.2f} ± {Hoehe23_err*1e9:.2f} nm')
print(f'Höhe 12nm: {Hoehe12*1e9:.2f} ± {Hoehe12_err*1e9:.2f} nm')
