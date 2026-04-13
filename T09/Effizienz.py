#%%

import numpy as np
import uncertainties as unc
from uncertainties import ufloat
from uncertainties import unumpy as unp
import math
N_H = 10000



def split_uncertainty(result, stat_tag="stat", syst_tag="syst"):
    """
    Zerlegt die Unsicherheit eines uncertainties-Ergebnisses
    in statistische und systematische Beiträge anhand der Tags.
    """
    components = result.error_components()
    stat2 = 0.0
    syst2 = 0.0
    for var, err in components.items():
        if var.tag is None:
            continue
        if stat_tag in var.tag:
            stat2 += err**2
        elif syst_tag in var.tag:
            syst2 += err**2
    return (
        result.nominal_value,
        math.sqrt(stat2),
        math.sqrt(syst2)
    )

N_data_H = np.array([
    [3910, 3945, 3938, 3973],
    [1670, 1685, 1716, 1731],
    [1867, 1927, 1966, 2028]
])

N_MC_H = np.array([
    [9791, 9863, 9832, 9905],
    [9452, 9536, 9607, 9693],
    [8719, 8972, 9174, 9431]
])


N_data_mu = np.array([
    [98-1, 106, 99, 107], #-1 wegen Fehlmessung
    [26, 31, 29, 35],
    [30, 39, 39, 48]
])
N_mu = 9969

N_MC_mu = np.array([
    [5272, 5452, 5310, 5400],
    [4593, 4830, 4801, 5050],
    [4615, 4920, 4731, 5050]
])


#systematische_Fehler = [21, 40.87, 43.21]
systematische_Fehler = [3, 2.83, 6]
eff_H =  N_MC_H[0, 0] / N_H

eff_mu = N_MC_mu / N_mu
labels = [91, 89, 93]
for i in range(3):
    print(f"statistische Fehler: N_nominal {labels[i]}: {np.sqrt(N_data_mu[i, 0]):.2f}")
    syst_Fehler = 0
    for j in range(4):
        syst_Fehler += np.abs(N_data_mu[i, j] - N_data_mu[i, 0])/2
    syst_Fehler /= 3
    print(f"systematische Fehler: N_nominal {labels[i]}: {syst_Fehler:.2f}")
print('---'*20)
#eff_n = np.array([.9791, .9452, 0.8719])        
eff_n = np.array([.5288, .4607, 0.4629])        
eff_std = np.sqrt(eff_n * (1 - eff_n) / N_H)
#print("Effizienz und Fehler:", eff_n*100, eff_std*100)
for i in range(3):
    N_stat = ufloat(N_data_mu[i, 0], np.sqrt(N_data_mu[i, 0]), tag=f"stat_N")
    N_syst = ufloat(0, systematische_Fehler[i], tag=f"syst_N")
    N = N_stat + N_syst
    eff = ufloat(eff_n[i], eff_std[i], tag=f"stat_eff")
    A = N / eff
    value, stat, syst = split_uncertainty(A)
    print(f"N_nominal {labels[i]}: A = {value:.3f} ± {stat:.3f} (stat) ± {syst:.3f} (syst)")    

L = np.array([134.4, 175.4, 151.1])

for i in range(3):
    N_stat = ufloat(N_data_mu[i, 0], np.sqrt(N_data_mu[i, 0]), tag=f"stat_N")
    N_syst = ufloat(0, systematische_Fehler[i], tag=f"syst_N")
    N = N_stat + N_syst
    eff = ufloat(eff_n[i], eff_std[i], tag=f"stat_eff")
    L_ufloat = ufloat(L[i], 0.01*L[i], tag="syst_L")
    A = N / eff
    sigma = A / L_ufloat
    value, stat, syst = split_uncertainty(sigma)
    print(f"N_nominal {labels[i]}: crossection = {value:.3f} ± {stat:.3f} (stat) ± {syst:.3f} (syst)")
    print(f"Gesamtfehler: {np.sqrt(stat**2 + syst**2):.3f}")
for (var, error) in sigma.error_components().items():
    print( "{}: {}".format(var.tag, error))

eeff_mu = N_MC_mu[2, 0] / N_mu
eeff_mu_std = np.sqrt(eeff_mu * (1 - eeff_mu) / N_mu)

print('----'*20)

Gamma_e = ufloat(83.91, 0.12)*1e-3  # MeV
Gamma_myons = ufloat(83.99, 0.18)*1e-3  # MeV
mZ = ufloat(91.1876, 0.0021)  # GeV
Gamma_Z = ufloat(2.4952, 0.0023)  # GeV
sigma_0 = 12 * np.pi * (Gamma_e * Gamma_myons) / (mZ**2 * Gamma_Z**2)*(389379)  # nb
print(sigma_0)
print('----'*20)
#mZ = ufloat(91.17, 0.01)
mZ = ufloat(91.17, 0.01)
Gamma_Z= ufloat(2.535, 0.047)
sigma_0 = ufloat(1.574, 0.145)
Gamma_e = unp.sqrt(sigma_0/(389379) * mZ**2 * Gamma_Z**2 / (12 * np.pi))
print(f"Gamma_e1: {Gamma_e*1000:.3f} MeV")
print('----'*20)
G_F = 1.166e-5  # GeV^-2
Gamma_v = G_F*mZ**3 / (24*np.sqrt(2) * np.pi) * 2
B =Gamma_Z/3-Gamma_v
sigma_0_had = ufloat(40.09, 0.72)/ (389379)
C = mZ**2 / (12 * np.pi) * Gamma_Z**2 *sigma_0_had /3
print(C, (B/2)**2)
Gamma_e = B/2-unp.sqrt((B/2)**2 - C)
print(f"Gamma_e2: {Gamma_e*1000:.3f} MeV")
print('----'*20)
A = Gamma_e*24*np.sqrt(2)*np.pi/(G_F*mZ**3)
sin_t_w = 1/4*(1+unp.sqrt(A+A.s-1))
#print(f"sin^2(theta_W): {sin_t_w:.3f}")
sin_t_w = 1/4*(1-unp.sqrt(A+A.s-1))

print(f"sin^2(theta_W): {sin_t_w:.3f}")
sin_t_w2 = 1/4*(1-unp.sqrt(A+2*A.s-1))
print(f"unterer Fehler: {sin_t_w2-sin_t_w:.3f}")

sin_t_w = ufloat(0.232, 0.025)
gamma_d = G_F*mZ**3 / (24*np.sqrt(2) * np.pi) * (1+(1-4*1/3*sin_t_w)**2)
gamma_u = G_F*mZ**3 / (24*np.sqrt(2) * np.pi) * (1+(1-4*2/3*sin_t_w)**2)

gamma_had = 1.04*(2*gamma_u + 3*gamma_d)*1000
print(f"Gamma_had: {gamma_had:.3f} GeV")
gamma_lit = ufloat1744 = ufloat(1744.4, 2.0)
print(f"Gamma_had / Gamma_lit: {gamma_lit/gamma_had:.3f}")
# %%
