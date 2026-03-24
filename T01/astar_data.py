#%%

import numpy as np
from sympy import E
from uncertainties import ufloat
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


rho = 1.2041*1e-3

data = np.loadtxt('astar_data.txt', skiprows = 9)

def linear(x, m, b):
    return m * x + b
Kanalnummern = [7794, 7400, 6384, 900]
Kanalnummern_unc = [1, 25, 50, 50]
proj_range = data[:,2]/rho
stopping_power = data[:,1]*rho
restreichweite = np.array([7.06, 7.06-1, 7.06-1.8, 7.06-4.18])
restreichweite = restreichweite-2.95
restreichweite_unc = np.array([0.35, 0.35, 0.35, 0.35])
Es = []
for i in range(len(restreichweite)):
    R = ufloat(restreichweite[i], restreichweite_unc[i])
    E = np.argmin(np.abs(proj_range - R))
    E = data[E,0]
    print(f"Restreichweite: {R:.2f} cm, Energie: {E:.2f} MeV")
    Es.append(E)
popt, pcov = curve_fit(linear, Es, Kanalnummern, sigma=Kanalnummern_unc, absolute_sigma=True, p0 =[900, 1500])
plt.plot(Es, Kanalnummern, 'o', label='Datenpunkte')
plt.plot(Es, linear(np.array(Es), *popt), label='Lineare Anpassung')


Peakpos = [2354.2, 3549.1, 3994.8, 4891.3, 7793.7]
peakpos_unc = [15, 33, 21, 6, 1]

Reichweite_lit = np.array([3.36, 3.9, 3.93, 4.15, 7.06]) 
Reichweite_lit = Reichweite_lit-2.95

Es = []
for i in range(len(Reichweite_lit)):
    R = Reichweite_lit[i]
    E = np.argmin(np.abs(proj_range - R))
    E = data[E,0]
    print(f"Restreichweite: {R:.2f} cm, Energie: {E:.2f} MeV")
    Es.append(E)
popt, pcov = curve_fit(linear, Es, Peakpos, sigma=peakpos_unc, absolute_sigma=True, p0 =[900, 1500])
plt.plot(Es, Peakpos, 'o', label='Datenpunkte')
plt.plot(Es, linear(np.array(Es), *popt), label='2')
plt.legend()