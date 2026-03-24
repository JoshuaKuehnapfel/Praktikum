#%%

import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit
from uncertainties import ufloat
intervall = np.array([1000, 1000, 1000, 1000])
def multi_gaussian(x, A1, A2, A3, A4, x1, x2, x3, x4, sigma1, sigma2, sigma3, sigma4):
    return (
        A1 * np.exp(-0.5 * ((x - x1) / sigma1) ** 2) +
        A2 * np.exp(-0.5 * ((x - x2) / sigma2) ** 2) +
        A3 * np.exp(-0.5 * ((x - x3) / sigma3) ** 2) +
        (A4 * np.exp(-0.5 * ((x - x4) / sigma4) **2 if A4 != 0 and x4 != 0 else 0))
    )
def single_gaussian(x, A, x0, sigma):
    return A * np.exp(-0.5 * ((x - x0) / sigma) ** 2)
N_lower = 800

def peak_fitter(data, x_data, p0_array, axis):
    popt, pcov = curve_fit(multi_gaussian, x_data, data, p0=p0_array, sigma = np.sqrt(data+0.001), absolute_sigma=True, maxfev=10000)
    x_smooth = np.linspace(min(x_data), max(x_data), 2000)
    axis.plot(x_smooth, multi_gaussian(x_smooth, *popt), color='black', label='Fit', zorder = 2)
    for i in range(4):
        peak_pos = popt[4 + i]
        peak_amp = popt[i]
        peak_sigma = popt[8 + i]
        x_peak = np.linspace(peak_pos - 1500, peak_pos + 1500, 500)
        axis.plot(x_peak, single_gaussian(x_peak, peak_amp, peak_pos, peak_sigma), linestyle='--', label=f'Peak {i+1}', zorder=2)
    return popt, pcov


daten_namen = ['Spektrum','1.Peak_links','1.Peak_Mitte','1.Peak_rechts',
               '2.Peak_links', '2.Peak_Mitte', '2.Peak_rechts',
               '3.Peak_links', '3.Peak_Mitte', '3.Peak_rechts',
               '4.Peak_links', '4.Peak_Mitte', '4.Peak_rechts'] 

peaks ={
    'Spektrum': [1800, 3700, 5000, 7800],
    '1.Peak_links': [1250, 3500, 7600,0],
    '1.Peak_Mitte': [800, 3500, 7400,0],
    '1.Peak_rechts': [900, 2700, 7200,0],
    '2.Peak_links': [1250, 3500, 7600,0],
    '2.Peak_Mitte': [800, 3500, 7400,0],
    '2.Peak_rechts': [900, 2700, 7200,0],
    '3.Peak_links': [1250, 3500, 7600,0],
    '3.Peak_Mitte': [800, 3500, 7400,0],
    '3.Peak_rechts': [900, 2700, 7200,0],
    '4.Peak_links': [1250, 3500, 7600,0],
    '4.Peak_Mitte': [800, 3500, 7400,0],
    '4.Peak_rechts': [900, 2700, 7200,0]

}


daten = np.loadtxt('Alpha in Luft/Alpha_Luft_'+daten_namen[0]+'.TKA')
d_name = daten_namen[0]
d = daten[N_lower:9000]
daten_uns = np.sqrt(d)
x = np.arange(len(d))+N_lower
plt.clf()
fig, ax = plt.subplots(2,1,sharex=True, figsize=(8,6), gridspec_kw={'height_ratios':[4,2]})
ax[0].errorbar(x, d, yerr=daten_uns,fmt = '.', alpha=1, label='Daten', zorder = 1)
fig.suptitle(f"Gesamtspektrum mit Fit")
p0 = [peaks[d_name][0], peaks[d_name][1], peaks[d_name][2], peaks[d_name][3], peaks[d_name][0], peaks[d_name][1], peaks[d_name][2], peaks[d_name][3], 1, 1, 1, 1]
popt, pcov = peak_fitter(d, x, p0, ax[0])
ax[0].vlines(popt[4:8], ymin=0, ymax=max(d), color='red', label='Peakpositionen', zorder = 3, linestyle='dashdot')
print(f"Peak positionen: {popt[4:8]}+-{np.sqrt(np.diag(pcov))[4:8]}")
residuen = d - multi_gaussian(x, *popt)
chiq = np.sum(residuen**2 / (daten_uns**2+0.001))
ax[1].errorbar(x, residuen,yerr=daten_uns,fmt='.',linestyle='none',color='maroon',label=r'$\chi^2$ ' + f'= {chiq:.3f}')
print(f"Reduziertes $\chi^2$: {chiq/(len(d)-12):.3f}", f"Anzahl Datenpunkte: {len(d)}", f"Anzahl Fit-Parameter: 12")
ax[1].axhline(0, min(x), max(x), linewidth = 2, color ='black')
for a in ax:
    a.grid()
    a.legend()
#plt.savefig('Gesamtspektrum_mit_Fit.svg', dpi=300)



fig, ax =plt.subplots(2, 1,sharex=True, figsize=(8,6), dpi=300, gridspec_kw={'height_ratios':[4,2]})
Energien = np.array([4.69, 5.3, 6.00, 7.69])
x_smooth = np.linspace(min(Energien), max(Energien), 100)
def linear(x, m, b):
    return m * x + b
sigmas = np.sqrt(np.diag(pcov))[4:8]
ax[0].errorbar(Energien, popt[4:8], yerr = sigmas, fmt = 'o', label='Peakpositionen')
popt2, pcov = curve_fit(linear,Energien, popt[4:8], sigma=sigmas, absolute_sigma=True)
ax[0].plot(x_smooth, linear(x_smooth, *popt2), color='black', label='Fit', zorder=2)
ax[0].set_xlabel('Energie (MeV)')
ax[0].set_ylabel('Position (channels)')
ax[1].errorbar(Energien, popt[4:8] - linear(Energien, *popt2), yerr=sigmas, fmt='o', color='maroon', label='Residuen')
ax[1].axhline(0, min(Energien), max(Energien), linewidth=2, color='black')
ax[1].set_xlabel('Energie (MeV)')   
ax[1].set_ylabel('Residuen')
for a in ax:
    a.grid()
    a.legend()
fig.suptitle('Kalibrierung mit linearem Fit')
plt.savefig('Kalibrierung_mit_Fit.svg', dpi=300)

print(f"Kalibrierungsparameter: m = {popt2[0]:.3f} +- {np.sqrt(pcov[0,0]):.3f}, b = {popt2[1]:.3f} +- {np.sqrt(pcov[1,1]):.3f}")
print(f"Reduziertes $\chi^2$: {np.sum(((popt[4:8] - linear(Energien, *popt2)) / sigmas) ** 2) / (len(Energien) - 2):.3f}")
print(np.sum(((popt[4:8] - linear(Energien, *popt2)) / sigmas) ** 2))

Abstand = np.array([27.9, 27.3, 26.5])
Abstand = Abstand[0]-Abstand
Abstand_unc = [0.1, 0.1, 0.1]
fig, ax = plt.subplots(2, 1,sharex=True, figsize=(8,6), dpi=300, gridspec_kw={'height_ratios':[4,2]})
ax[1].plot(Abstand, Energien[:-1])
plt.show()

# %%
