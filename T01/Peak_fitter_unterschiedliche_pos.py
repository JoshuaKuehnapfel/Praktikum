#%%

import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit
from sympy import residue

N_lower = 800
intervall = np.array([1000, 1000, 1000, 1000])
def single_gaussian(x, A, x0, sigma):
    return A / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - x0) / sigma) ** 2)
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
        print(f"Fläche unter Kurve : {peak_amp:.2f} +- {np.sqrt(pcov[i,i]):.2f}")
    return popt, pcov

daten_namen = ['Spektrum','1.Peak_links','1.Peak_Mitte','1.Peak_rechts',
               '2.Peak_links', '2.Peak_Mitte', '2.Peak_rechts',
               '3.Peak_links', '3.Peak_Mitte', '3.Peak_rechts']
               #'4.Peak_links', '4.Peak_Mitte', '4.Peak_rechts'] 

peaks ={
    'Spektrum': [1800, 3700, 5000, 7800],
    '1.Peak_links': [1250, 3500, 7750],
    '1.Peak_Mitte': [800, 3500, 7750],
    '1.Peak_rechts': [900, 2700, 7650],
    '2.Peak_links': [1250, 3500, 7400],
    '2.Peak_Mitte': [800, 3500, 7300],
    '2.Peak_rechts': [900, 2700, 7400],
    '3.Peak_links': [1250, 3500, 6500],
    '3.Peak_Mitte': [800, 3500, 6500],
    '3.Peak_rechts': [900, 2700, 6300],
    '4.Peak_links': [1250, 3500, 1000],
    '4.Peak_Mitte': [800, 3500, 1000],
    '4.Peak_rechts': [900, 2700, 1000]

}
A = []
for d_name in daten_namen:
    daten = np.loadtxt('Alpha in Luft/Alpha_Luft_'+d_name+'.TKA')
    d = daten[N_lower:9000]
    daten_uns = np.sqrt(d)
    x = np.arange(len(d))+N_lower
    plt.clf()
    fig, ax = plt.subplots(2,1,sharex=True, figsize=(8,6), gridspec_kw={'height_ratios':[3,1], 'hspace':0.05}, layout = 'tight')
    ax[0].errorbar(x, d, yerr=daten_uns,fmt = '.', alpha=1, label='Daten', zorder = 1)
    p0 = [max(d), peaks[d_name][-1], 1]
    popt, pcov = curve_fit(single_gaussian, x[peaks[d_name][-1]-1500:peaks[d_name][-1]+1500], d[peaks[d_name][-1]-1500:peaks[d_name][-1]+1500], p0=p0, sigma = np.sqrt(d+0.001)[peaks[d_name][-1]-1500:peaks[d_name][-1]+1500], absolute_sigma=True, maxfev=10000)
    x_smooth = np.linspace(popt[1] - 1500, popt[1] + 1500, 2000)
    mask = (x >= x_smooth[0]) & (x <= x_smooth[-1])
    x_res = x[mask]
    d_res = d[mask]
    print(f"Fläche unter Kurve : {popt[0]:.2f} +- {np.sqrt(pcov[0,0]):.2f}")
    A.append(popt[0])
    daten_uns_res = daten_uns[mask]
    residuen = d_res - single_gaussian(x_res, *popt)
    ax[0].plot(x_smooth, single_gaussian(x_smooth, *popt), color='black', label='Fit', zorder = 2)
    ax[0].vlines(popt[1], min(d), max(d), color='maroon', linestyle='--', label=f'Peakposition: {popt[1]:.0f}'r' $\pm$ 'f'{np.sqrt(np.diag(pcov))[1]:.0f}' , zorder=2)  
    ax[0].set_title(f"{d_name} mit Fit")
    ax[0].set_ylabel('Count')
    ax[0].grid(True, which='both', linestyle='--', linewidth=0.5)
    ax[0].legend()
    ax[1].set_xlabel('Kanäle')
    ax[1].set_ylabel('Count')
    ax[1].vlines(popt[1], min(residuen), max(residuen), color='maroon', linestyle='--',  zorder=2)  
    ax[1].grid(True, which='both', linestyle='--', linewidth=0.5)
    ax[1].errorbar(x_res, residuen, yerr=daten_uns_res, fmt='.', color='lightskyblue', label=r'Residuen', zorder=1)
    ax[1].axhline(0, linewidth = 1.5, linestyle = 'dotted', color ='black', zorder = 3)
    chisquare = np.sum((residuen)**2 / (daten_uns_res+0.1**2))
    dof = len(residuen) - len(popt)
    ax[1].text(0.05, 0.95, r'$\frac{\chi^2}{Ndof}$ ' + f'= {chisquare:.2f}' + f', Ndof={dof}', transform=ax[1].transAxes, verticalalignment='top')
    ax[1].text(0.05, 0.75, r'$ \Rightarrow \quad \frac{\chi^2}{Ndof}$ ' + f'= {chisquare/dof:.2f}', transform=ax[1].transAxes, verticalalignment='top')
    ax[1].legend()
    #plt.savefig(f'Alpha in Luft/Plots/{d_name}_fit.png', dpi=300)
    #print(f"{d_name} - Peakposition: {popt[1]:.0f} +- {np.sqrt(np.diag(pcov))[1]:.0f}, Reduziertes chi^2: " + r'$\frac{\chi^2}{N_{dof}}$ ' + f'= {chisquare/dof:.2f}')
    #plt.show()
Abstand = np.array([28.29,27.9, 27.3, 26.5, 24.1])
Abstand = Abstand[0]-Abstand
Abstand_unc = [0.05,0.05, 0.05, 0.05, 0.1]
plt.clf()
def linear(x, m, b):
    return m * x + b
Abstand = Abstand[1:]
Abstand_unc = Abstand_unc[1:]

peak_pos = [7668, 7244, 6384, 900]
peak_pos_unc = [25, 40, 50, 50]
peak_pos_unc_ges = np.sqrt(np.array(peak_pos_unc)**2) 
peaks_ges = peak_pos

sigma2 = np.sqrt(np.array(peak_pos_unc_ges)**2 + np.array(Abstand_unc)**2)
popt, pcov = curve_fit(linear, Abstand, peaks_ges, sigma=sigma2, absolute_sigma=True)
x_smooth = np.linspace(min(Abstand), max(Abstand)+0.1, 100)
fig, ax = plt.subplots(2, 1, figsize=(8,6), gridspec_kw={'height_ratios': (3,1), 'hspace': 0.05}, dpi=300, sharex=True)
ax[0].plot(x_smooth, linear(x_smooth, *popt), color='maroon', label='Fit', zorder=2)
ax[0].set_ylabel('Peakposition (channels)')
ax[0].errorbar(Abstand, peaks_ges, xerr=Abstand_unc, yerr=peak_pos_unc_ges, fmt='o', color='black', capsize=5, label='Peakpositionen', zorder=2)
ax[1].set_ylabel('Residuen')
residuen = peaks_ges - linear(Abstand, *popt)
ax[1].errorbar(Abstand, residuen, xerr=Abstand_unc, yerr=peak_pos_unc_ges, fmt='o', color='maroon', capsize=5, label='Residuen', zorder=2)
ax[1].axhline(0, linewidth=1.5, linestyle='--', color='black', zorder=3)
chisquare = np.sum((residuen)**2 / (sigma2**2))
ax[0].text(0.05, 0.5, f'm = {popt[0]:.1f}' + r' $\pm$ ' + f'{np.sqrt(pcov[0,0]):.1f}' + '\n' + f'b = {popt[1]:.1f}' + r' $\pm$ ' + f'{np.sqrt(pcov[1,1]):.1f}', transform=ax[0].transAxes, verticalalignment='top', fontsize=12,bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))  
ax[1].text(0.25, 0.45, r'$\frac{\chi^2}{N_{dof}}$ ' + f'= {chisquare/2:.2f}' + f', Ndof={len(peaks_ges)-len(popt)}', transform=ax[1].transAxes, verticalalignment='top', fontsize=12,bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
ax[1].set_xlabel(r'$\Delta$ Abstand (cm)')
for a in ax:
    a.grid(True, which='both', linestyle='--', linewidth=0.5)
    a.legend()
fig.suptitle('Peakpositionen des 7.69 MeV-Peaks in Abhängigkeit des Abstandes')
plt.savefig('Alpha in Luft/Plots/Peakpositionen_Abstand.svg', dpi=300)

Rauschen = np.loadtxt('Alpha in Luft/Alpha_Luft_Rausch_4,4.TKA')

peak_4_names = ['4.Peak_links', '4.Peak_Mitte', '4.Peak_rechts']
for name in peak_4_names:
    daten = np.loadtxt('Alpha in Luft/Alpha_Luft_'+name+'.TKA')-Rauschen
    d = daten[N_lower:1800]
    A.append(d.sum())
Abstand = np.array([28.29,27.9,27.8, 27.7, 27.4, 27.3,27.2,26.6, 26.5,26.4, 24.2, 24.1, 24.0])
Messzeiten = np.array([300, 90, 90, 90, 90, 90, 90,90, 90, 90, 120, 120, 120])
A = np.array(A)
Abstand = Abstand[0]-Abstand
Abstand_unc = np.ones_like(Abstand)*0.01
Abstand_unc[-3:] = 0.05
Raten = A / Messzeiten
Raten_unc = np.sqrt(A) / Messzeiten


plt.clf()

popt, pcov = curve_fit(linear, Abstand[-3:], Raten[-3:], sigma=Raten_unc[-3:], absolute_sigma=True)
x_smooth = np.linspace(max(Abstand)-0.8, max(Abstand)+0.1, 100)
max_reichweite = (0 - popt[1]) / popt[0]
max_reichweite_unc = max_reichweite * np.sqrt((np.sqrt(pcov[0,0])/popt[0])**2 + (np.sqrt(pcov[1,1])/popt[1])**2)
plt.plot(x_smooth, linear(x_smooth, *popt), color='Blue', label='Fit', zorder=2)
plt.axvline(max_reichweite, color='grey', linestyle='--', label='maximale Reichweite', zorder=2)
plt.axvspan(max_reichweite - max_reichweite_unc, max_reichweite + max_reichweite_unc, alpha=0.2, color='grey', zorder=1)
plt.errorbar(Abstand, Raten, xerr=Abstand_unc, yerr=Raten_unc, fmt='o', color='maroon', capsize=5, label='Raten', zorder=2)
plt.xlabel(r'$\Delta$ Abstand (cm)')
plt.ylabel('Rate (counts/s)')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.title('Raten in Abhängigkeit des Abstandes')
plt.legend()
plt.text(0.3, 0.7, 'Die maximale Reichweite der 7.69 MeV\n' r'$\alpha$-Teilchen liegt bei ' + f'{max_reichweite:.2f} $\pm$ {max_reichweite_unc:.2f} cm', transform=plt.gca().transAxes, verticalalignment='top', fontsize=12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
plt.savefig('Alpha in Luft/Plots/Raten_Abstand.svg', dpi=300)



# %%
