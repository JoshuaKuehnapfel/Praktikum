#%%

import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy as unp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

abstand = np.array([
    9.05, 8.9, 8.7, 8.5, 8.3, 8.1, 7.9,
    7.7, 7.5, 7.0, 6.5, 6.0, 5.0, 4.5
])

strom = np.array([
    2.90, 2.68, 2.68, 2.40, 2.04, 1.68, 1.32,
    1.02, 0.775, 0.435, 0.365, 0.360, 0.115, 0.0125
])

strom_err = np.array([
    0.03, 0.03, 0.03, 0.05, 0.04, 0.20, 0.20,
    0.15, 0.015, 0.015, 0.015, 0.015, 0.010, 0.0015
])
system_fehler = 0.03*strom
strom_err = np.sqrt(strom_err**2 + system_fehler**2)

abstand = abstand[0] - abstand
strom_unc = unp.uarray(strom, strom_err)


def I_x(x, I_0, alpha):
    return I_0 * np.exp(-alpha * x)


def plot_errorfit(x, y, y_err, c='orange', p0 = (1, 1), axis = plt.gca()):
    popt, pcov = curve_fit(I_x, x, y, sigma=y_err, absolute_sigma=True, p0 = p0)
    x_fit = np.linspace(min(x), max(x), 300)
    y_fit = I_x(x_fit, *popt)
    I0, alpha = popt
    sigma = np.sqrt(
        np.exp(2*-alpha*x_fit)*pcov[0,0]
        + (I0*x_fit*np.exp(-alpha*x_fit))**2 * pcov[1,1]
        - 2*np.exp(-alpha*x_fit)*(I0*x_fit*np.exp(-alpha*x_fit))*pcov[0,1]
    )
    axis.plot(x_fit, y_fit)
    axis.fill_between(x_fit, y_fit-3*sigma, y_fit+3*sigma, alpha=0.3, color=c, label=r'3$\sigma$  Unsicherheitsbereich')
    axis.plot(x_fit, y_fit, label='Exponentielle Anpassung', color=c)
    axis.set_title('Strom als Funktion des Abstands')
    residuen = (y - I_x(x, *popt))
    print(f'Angepasste Parameter: I0 = {I0:.3f} nA +- {np.sqrt(pcov[0,0]):.3f}, alpha = {alpha:.3f} +- {np.sqrt(pcov[1,1]):.3f} 1/cm')
    return residuen
fig, ax = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(8, 6.5))
Gl = 3
GU = -4
ax[0].errorbar(abstand, unp.nominal_values(strom_unc), yerr=unp.std_devs(strom_unc), fmt='o', label='Daten mit Fehlern')
residuen = plot_errorfit(abstand[Gl:GU], unp.nominal_values(strom_unc[Gl:GU]), unp.std_devs(strom_unc[Gl:GU]), c='orange', p0=(3, 0.5), axis=ax[0])
chisq = np.sum((residuen / unp.std_devs(strom_unc[Gl:GU]))**2)

plateu_idx = [(1,2), (10, 11)]
plateu_labels = ['Plateau 1', 'Plateau 2']
plateau_unc_labels = [r'1 $\sigma$ Unsicherheitsbereich', '']
plateu_c = ['limegreen', 'darkgreen']
for idx in plateu_idx:
    j = plateu_idx.index(idx)
    höhe = np.average(unp.nominal_values([strom_unc[i]for i in idx]), weights=1/unp.std_devs([strom_unc[i]for i in idx])**2)
    höhe_err = np.sqrt(1/np.sum(1/unp.std_devs([strom_unc[i]for i in idx])**2))
    ax[0].axhline(höhe, color=plateu_c[j], linestyle='--', label=plateu_labels[j])
    ax[0].fill_between(abstand, höhe-höhe_err, höhe+höhe_err, alpha=0.3, color=plateu_c[j], label=plateau_unc_labels[j])
    ax[1].errorbar(abstand[idx[0]:idx[-1]+1], unp.nominal_values([strom_unc[i] for i in idx]) - höhe, yerr=unp.std_devs([strom_unc[i] for i in idx]), fmt='o', color=plateu_c[j], label=f'Residuen {plateu_labels[j]}')

ax[0].set_ylabel('Strom [nA]')
ax[0].legend() 
ax[1].errorbar(abstand[Gl:GU], residuen, yerr=unp.std_devs(strom_unc[Gl:GU]), fmt='o', color='red', linestyle='None', label=r'Residuen ($\chi^2$ = {:.2f})'.format(chisq))
ax[1].set_xlabel('Abstand [cm]')
ax[1].legend()
ax[1].set_ylabel('Residuen [nA]')
ax[1].axhline(0, color='black', linestyle='--')

for a in ax:
    a.grid()
fig.tight_layout()
plt.savefig('Bilder/Absorptionskoeffizient.svg', dpi=300)
# %%
