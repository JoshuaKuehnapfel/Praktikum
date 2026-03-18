#%%

import numpy as np
from sympy import E
from uncertainties import ufloat
from scipy.optimize import curve_fit

rho = 1.2041*1e-3

data = np.loadtxt('astar_data.txt', skiprows = 9)

proj_range = data[:,2]/rho
stopping_power = data[:,1]*rho
restreichweite = np.array([7.06, 6.08, 5.27, 2.87])
restreichweite_unc = np.array([0.35, 0.35, 0.35, 0.35])
for i in range(len(restreichweite)):
    R = ufloat(restreichweite[i], restreichweite_unc[i])
    E = np.argmin(np.abs(proj_range - R))
    print(E)
    E = data[E,0]
    print(f"Restreichweite: {R:.2f} cm, Energie: {E:.2f} MeV")
