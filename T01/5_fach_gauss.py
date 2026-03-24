#%%
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit
from uncertainties import ufloat

N_lower = 800

def multi_gaussian(x, A1, A2, A3, A4, A5, x1, x2, x3, x4, x5, sigma1, sigma2, sigma3, sigma4, sigma5):
    return (
        A1 * np.exp(-0.5 * ((x - x1) / sigma1) ** 2) +
        A2 * np.exp(-0.5 * ((x - x2) / sigma2) ** 2) +
        A3 * np.exp(-0.5 * ((x - x3) / sigma3) ** 2) +
        (A4 * np.exp(-0.5 * ((x - x4) / sigma4) **2 if A4 != 0 and x4 != 0 else 0) +
         A5 * np.exp(-0.5 * ((x - x5) / sigma5) **2 if A5 != 0 and x5 != 0 else 0))
    )
def single_gaussian(x, A, x0, sigma):
    return A * np.exp(-0.5 * ((x - x0) / sigma) ** 2)


def peak_fitter(data, x_data, p0_array, axis):
    popt, pcov = curve_fit(multi_gaussian, x_data, data, p0=p0_array, sigma = np.sqrt(data+0.001), absolute_sigma=True, maxfev=10000, bounds=([0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 1, 1, 1, 1, 1], [np.inf, np.inf, np.inf, np.inf, np.inf, 9000, 9000, 9000, 9000, 9000, np.inf, np.inf, np.inf, np.inf, np.inf]))
    x_smooth = np.linspace(min(x_data), max(x_data), 2000)
    axis.plot(x_smooth, multi_gaussian(x_smooth, *popt), color='black', label='Fit', zorder = 2)
    for i in range(5):
        peak_pos = popt[5 + i]
        peak_amp = popt[i]
        peak_sigma = popt[10 + i]
        x_peak = np.linspace(peak_pos - 1500, peak_pos + 1500, 500)
        axis.plot(x_peak, single_gaussian(x_peak, peak_amp, peak_pos, peak_sigma), linestyle='--', label=f'Peak {i+1}', zorder=2)
    return popt, pcov


daten_namen = ['Spektrum','1.Peak_links','1.Peak_Mitte','1.Peak_rechts',
               '2.Peak_links', '2.Peak_Mitte', '2.Peak_rechts',
               '3.Peak_links', '3.Peak_Mitte', '3.Peak_rechts',
               '4.Peak_links', '4.Peak_Mitte', '4.Peak_rechts'] 

peaks ={
    'Spektrum': [1800, 3700, 4100, 5000, 7800],
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
fig.suptitle(f"Gesamtspektrum mit Fit (5 Peaks)")
p0 = [peaks[d_name][0], peaks[d_name][1], peaks[d_name][2], peaks[d_name][3], peaks[d_name][4], peaks[d_name][0], peaks[d_name][1], peaks[d_name][2], peaks[d_name][3], peaks[d_name][4], 1, 1, 1, 1, 1]
popt, pcov = peak_fitter(d, x, p0, ax[0])
ax[0].vlines(popt[5:10], ymin=0, ymax=max(d), color='red', label='Peakpositionen', zorder = 3, linestyle='dashdot')
print(f"Peak positionen: {popt[5:10]}+-{np.sqrt(np.diag(pcov))[5:10]}")
residuen = d - multi_gaussian(x, *popt)
chiq = np.sum(residuen**2 / (daten_uns**2+0.1**2))
print(len(d))
ax[1].errorbar(x, residuen,yerr=daten_uns,fmt='.',linestyle='none',color='maroon',label=r'$\chi^2$ ' + f'= {chiq:.3f}')
print(f"Reduziertes $\chi^2$: {chiq/(len(d)-15):.3f}", f"Anzahl Datenpunkte: {len(d)}", f"Anzahl Fit-Parameter: 15")
ax[1].axhline(0, min(x), max(x), linewidth = 2, color ='black')
for a in ax:
    a.grid()
    a.legend()
plt.savefig('Gesamtspektrum_mit_Fit_5.svg', dpi=300)





rho = 1.2041*1e-3
data = np.loadtxt('astar_data.txt', skiprows = 9)

proj_range = data[:,2]/rho
stopping_power = data[:,1]*rho

 
Peakpos = np.array([2354.2, 3549.1, 3994.8, 4891.3, 7793.7])
peakpos_unc = np.array([15.4, 33.7, 21.0, 6.0, 0.4])
Reichweite_lit = np.array([3.36, 3.9, 3.93, 4.15, 7.06]) 
Reichweite_lit = Reichweite_lit-2.95

E=[]
for i in range(len(Reichweite_lit)):
    idx = np.argmin(np.abs(proj_range - Reichweite_lit[i]))
    E.append(data[idx,0])
Energien = np.array(E)
Energien_unc = np.array([0.5, 0.5, 0.5, 0.5, 0.5])


fig, ax =plt.subplots(2, 1,sharex=True, figsize=(8,6), dpi=300, gridspec_kw={'height_ratios':[3,1], 'hspace':0.05}, layout = 'tight' )
x_smooth = np.linspace(min(Energien), max(Energien), 100)
def linear(x, m, b):
    return m * x + b
sigmas = np.sqrt(np.diag(pcov))[5:10]
ax[0].errorbar(Energien, popt[5:10], xerr = Energien_unc, yerr = sigmas, fmt = 'o', label='Peakpositionen')
popt2, pcov = curve_fit(linear,Energien, popt[5:10], sigma=sigmas, absolute_sigma=True)
ax[0].plot(x_smooth, linear(x_smooth, *popt2), color='black', label='Fit', zorder=2)
ax[0].set_xlabel(r'Energie $E_{detektiert}$ (MeV)')
ax[0].set_ylabel('Position (channels)')
ax[0].text(0.05, 0.95, f'm = {popt2[0]:.1f} +- {np.sqrt(pcov[0,0]):.1f}\nb = {popt2[1]:.1f} +- {np.sqrt(pcov[1,1]):.1f}', transform=ax[0].transAxes, verticalalignment='top', fontsize = 12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
ax[1].errorbar(Energien, popt[5:10] - linear(Energien, *popt2),xerr = Energien_unc,  yerr=sigmas, fmt='o', color='maroon', label='Residuen')
ax[1].axhline(0, linewidth=1.5, linestyle = 'dotted',color='black')
ax[1].set_xlabel(r'Energie $E_{detektiert}$ (MeV)')
ax[1].set_ylabel('Residuen')
ax[1].text(0.45, 0.95, f'Reduziertes $\chi^2$: {np.sum(((popt[5:10] - linear(Energien, *popt2))**2 / ((Energien_unc)**2+sigmas**2))) / (len(Energien) - 2):.3f}', transform=ax[1].transAxes, verticalalignment='top', fontsize=12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
for a in ax:
    a.grid()
    a.legend()
fig.suptitle('Kalibrierung mit linearem Fit')
plt.savefig('Kalibrierung_mit_Fit_5.svg', dpi=300)




# %%

