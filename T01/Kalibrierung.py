#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData
from scipy.optimize import curve_fit






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

def linear_func(p, x):
    m, c = p
    return m * x + c


linear_model = Model(linear_func)

datasc = RealData(Energien, Peakpos, sx=Energien_unc, sy=peakpos_unc)
odr = ODR(datasc, linear_model, beta0=[1000., 3000.])
output = odr.run()
output.pprint()
popt = output.beta
pcov = output.cov_beta






fig, ax =plt.subplots(2, 1,sharex=True, figsize=(8,6), dpi=300, gridspec_kw={'height_ratios':[3,1], 'hspace':0.05}, layout = 'tight' )
x_smooth = np.linspace(min(Energien), max(Energien), 100)
def linear(x, m, b):
    return m * x + b
popt2, pcov2 = curve_fit(linear, Energien, Peakpos, sigma=peakpos_unc, absolute_sigma=True)
ax[0].errorbar(Energien, Peakpos, xerr = Energien_unc, yerr = peakpos_unc, fmt = 'o', label='Peakpositionen')
ax[0].plot(x_smooth, linear(x_smooth, *popt2), color='black', label='Fit', zorder=2)
ax[0].set_xlabel(r'Energie $E_{detektiert}$ (MeV)')
ax[0].set_ylabel('Position (channels)')
ax[0].text(0.05, 0.95, f'm = {popt2[0]:.1f} +- {np.sqrt(pcov2[0,0]):.1f}\nb = {popt2[1]:.1f} +- {np.sqrt(pcov2[1,1]):.1f}', transform=ax[0].transAxes, verticalalignment='top', fontsize = 12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
residuals = Peakpos - linear(Energien, *popt2)
ax[1].errorbar(Energien, residuals, xerr = Energien_unc, yerr = peakpos_unc, fmt='o', color='maroon', label='Residuen')
ax[1].axhline(0, linewidth=1.5, linestyle = 'dotted',color='black')
ax[1].set_xlabel(r'Energie $E_{detektiert}$ (MeV)')
ax[1].set_ylabel('Residuen')
chi2_red = np.sum((residuals**2 / (peakpos_unc**2+Energien_unc**2))) / (len(Energien) - 2)
ax[1].text(0.45, 0.95, f'Reduziertes $\chi^2$: {chi2_red:.3f}', transform=ax[1].transAxes, verticalalignment='top', fontsize=12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
for a in ax:
    a.grid()
    a.legend()
fig.suptitle('Kalibrierung mit linearem Fit anhand des gesamten Spektrums')
plt.savefig('Kalibrierung_mit_Fit_5.svg', dpi=300)

#####


Kanalnummern = np.array([7794, 7400, 6384, 1000])
Kanalnummern_unc = np.array([1, 25, 50, 50])
proj_range = data[:,2]/rho
stopping_power = data[:,1]*rho
restreichweite = np.array([7.06, 7.06-1, 7.06-1.8, 7.06-4.18])
restreichweite = restreichweite-2.95
restreichweite_unc = np.array([0.35, 0.35, 0.35, 0.35])
Es = []
for i in range(len(restreichweite)):
    R = restreichweite[i]
    E = np.argmin(np.abs(proj_range - R))
    E = data[E,0]
    Es.append(E)
Energien = np.array(Es)
Energien_unc = np.array([0.5, 0.5, 0.5, 0.5])

datasc = RealData(Energien, Kanalnummern, sx=Energien_unc, sy=Kanalnummern_unc)
odr = ODR(datasc, linear_model, beta0=[1000., 3000.])
output = odr.run()
output.pprint()
popt = output.beta
pcov = output.cov_beta
popt, pcov = curve_fit(linear, Energien, Kanalnummern, sigma=Kanalnummern_unc, absolute_sigma=True)
x_smooth = np.linspace(min(Energien), max(Energien)+0.1, 100)
fig, ax = plt.subplots(2, 1, figsize=(8,6), gridspec_kw={'height_ratios': (3,1), 'hspace': 0.05}, dpi=300, sharex=True)
ax[0].plot(x_smooth, linear(x_smooth, *popt), color='maroon', label='Fit', zorder=2)
ax[0].set_ylabel('Peakposition (channels)')
ax[0].errorbar(Energien, Kanalnummern, xerr=Energien_unc, yerr=Kanalnummern_unc, fmt='o', color='black', capsize=5, label='Peakpositionen', zorder=2)
ax[1].set_ylabel('Residuen')
residuen = Kanalnummern - linear(Energien, *popt)
ax[1].errorbar(Energien, residuen, xerr=Energien_unc, yerr=Kanalnummern_unc, fmt='o', color='maroon', capsize=5, label='Residuen', zorder=2)
ax[1].axhline(0, linewidth=1.5, linestyle='--', color='black', zorder=3)
chisquare = np.sum((residuen)**2 / (Kanalnummern_unc**2))
ax[0].text(0.05, 0.5, f'm = {popt[0]:.1f}' + r' $\pm$ ' + f'{np.sqrt(pcov[0,0]):.1f}' + '\n' + f'b = {popt[1]:.1f}' + r' $\pm$ ' + f'{np.sqrt(pcov[1,1]):.1f}', transform=ax[0].transAxes, verticalalignment='top', fontsize=12,bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))  
ax[1].text(0.25, 0.45, r'$\frac{\chi^2}{N_{dof}}$ ' + f'= {chisquare/2:.2f}' + f', Ndof={len(Kanalnummern)-len(popt)}', transform=ax[1].transAxes, verticalalignment='top', fontsize=12,bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
ax[1].set_xlabel(r'Energie $E_{detektiert}$ (MeV)')
for a in ax:
    a.grid(True, which='both', linestyle='--', linewidth=0.5)
    a.legend()
fig.suptitle(r'Kalibrierung mit linearem Fit anhand der $^{214}$Po Peaks')
plt.savefig('Peakpositionen_Energie.svg', dpi=300)


# %%

