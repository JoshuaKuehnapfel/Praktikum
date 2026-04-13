#%%
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
Temperatur = np.array([15.1, 16, 18.3, 19.2, 20, 20.6, 21.2, 22.2, 22.9, 23.9, 24.8, 25.9, 27.3, 27.9, 28.9, 29.8, 30.7, 31.7, 32.3, 33, 34])
V = np.array([417*5, 4.10*5, 2.22*5, 2.1*5, 0.62*5, -9.5, -11, -11.1, -11.6, -12.1, -12.4, -12.4, -12.5, -12.8, -13.9, -14.3, -13.8, -13.1, -12.3, -10.4, -9.9])/50
V_uns = np.ones_like(V)*0.1/50
idx_cut = 5
plt.errorbar(Temperatur[idx_cut:], -V[idx_cut:], yerr=V_uns[idx_cut:], fmt='o')

plt.clf()
Tmin = 29.9
I = 600
T = Tmin - 0.004*(I-400)
print(f"Temperatur: {T:.2f} °C")
I = np.array([150,200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650])
Temperatur2 = np.array([30.9, 30.6, 30.7,30.5,  30.3, 30.1, 29.9, 29.7, 29.5, 29.3, 29.1, 28.9])
P = np.array([0, 0, 4, 8, 17, 25, 33, 41, 49, 57, 66, 76])*0.001 #in microwatt/cm**2
plt.plot(I, P, 'o', linestyle='dashed')
plt.clf()

I2 = np.array([200,225, 250, 275, 300, 350, 400, 450, 500, 550, 600])
P2 = np.array([0, 0, 0.22,0.95, 2.2, 4.67, 6.94, 9.22, 11.77, 14.07, 16.35])#in mW
P2unc = np.array([0, 0, 0.01, 0.01, 0.01, 0.01, 0.02, 0.01, 0.03, 0.04, 0.05])
plt.clf()
plt.plot(I2, P2, 'o', linestyle='dashed')
def linear(x, m, b):
    return m * x + b
def quadratic(x, a, b, c):
    return a * x**2 + b * x + c
def fit_mit_plot(func, xdata, ydata, sigma):
    popt, pcov = curve_fit(func, xdata, ydata, sigma=sigma, absolute_sigma=True)
    x_smooth = np.linspace(min(xdata), max(xdata), 100)
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8,6), dpi=300, gridspec_kw={'height_ratios': (3,1), 'hspace': 0.05})
    ax[0].errorbar(xdata, ydata, yerr=sigma, fmt='o', label='Leistung $P$ in mW', zorder=2)
    ax[0].plot(x_smooth, func(x_smooth, *popt), color='black', label='Fit', zorder=2)
    ax[0].set_xlabel('Strom $I$ in mA') 
    ax[0].set_ylabel('Leistung $P$ in mW')
    residuals = ydata - func(xdata, *popt)
    ax[1].errorbar(xdata, residuals, yerr=sigma, fmt='o', color='maroon', label='Residuen', zorder=2)
    ax[1].axhline(0, linewidth=1.5, linestyle='--', color='black', zorder=3)
    ax[1].set_xlabel('Strom $I$ in mA')
    ax[1].set_ylabel('Residuen')
    chi2 = np.sum((residuals / sigma) ** 2)
    ax[1].text(0.3, 0.45, r'$\frac{\chi^2}{N_{dof}}$ ' + f'= {chi2/(len(ydata)-len(popt)):.2f}' + f', Ndof={len(ydata)-len(popt)}', transform=ax[1].transAxes, verticalalignment='top', fontsize=12,bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    for a in ax:
        a.grid(True, which='both', linestyle='--', linewidth=0.5)
        a.legend()

fit_mit_plot(linear, I2, P2, P2unc)

I3 = np.array([625, 600, 550, 500, 450, 400, 350, 300])
P = np.array([0.354, 0.338, 0.250, 0.165, 0.106, 0.065, 0.028, 0.0057])#in mW
P_unc = np.array([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.00002])
plt.clf()
fit_mit_plot(quadratic, I3, P, P_unc)

Frequenz = np.array([429.60, 805.40, 916.0, 994.7, 1135.2, 1241.6, 1359.6, 1516.1, 1823.2, 1987.7, 2159.1, 2258.9,2458.7, 2588.0, 2697.8, 2831.8, 3008.7, 3143.5, 3418.0, 3861.2, 4414.4, 4214.5, 6366.6, 7930.5, 9137.3, 11205, 13434, 15192, 17237, 19575, 20206, 20502, 20632, 21047,21533])
Frequenz_unc = np.array([0.1, 0.05, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0.2, 0.2, 0.2,0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.7, 0.3, 0.3, 0.3, 1, 1, 1, 2, 2, 2, 2, 1,1, 1])
P3 = np.array([18, 36, 41, 45, 50, 54, 58, 65, 73, 77, 82, 88, 87, 85, 94, 91, 90, 92, 93, 97, 103, 99, 93, 102, 108, 94, 82, 86, 66, 74, 73, 58, 0, 0, 0])#in mikroW
P3_unc = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,5, 2, 3, 2, 3, 2, 2, 1, 1, 3, 2, 1, 1, 2, 1, 5, 3, 1, 2, 2,1, 0.1, 0.1,0.1])#in mikroW
plt.clf()
idx_cutoff = 12
fit_mit_plot(linear, Frequenz[:idx_cutoff], P3[:idx_cutoff], P3_unc[:idx_cutoff])
ax = plt.gcf().get_axes()
ax[0].plot(Frequenz[idx_cutoff:], P3[idx_cutoff:], 'o', linestyle='none', color = 'tab:blue')
ax[0].vlines(Frequenz[idx_cutoff], 0, P3.max(), colors='black', linestyles='dashed', alpha = 0.5)
ax[0].vlines(Frequenz[-3], 0, P3.max(), colors='black', linestyles='dashed', label='Grenzübergänge', alpha = 0.5)
ax[0].fill_betweenx([0, P3.max()+5], 0, Frequenz[idx_cutoff], color='blue', alpha=0.2, label='Lineare Region')
ax[0].fill_betweenx([0, P3.max()+5], Frequenz[idx_cutoff], Frequenz[-3], color='green', alpha=0.2, label='Nicht-lineare Region')
ax[0].fill_betweenx([0, P3.max()+5], Frequenz[-3], max(Frequenz)+800, color='orange', alpha=0.2, label='Sättigungsregion')
ax[0].legend(loc = 'lower center')
ax[0].set_xlabel('Frequenz [Hz]')
ax[0].set_ylabel(r'Leistung [$\mu$W]')
plt.savefig('Leistung_Frequenz.svg', dpi=300)
plt.clf()

Peak_shape = np.loadtxt('ALL0010/F0010CH1.csv', delimiter=',', skiprows=18, usecols=(3, 4))

def peak_fit(x, x0, l, A):
    return np.where(x > x0, A * np.exp(-l*(x-x0)), 0)

# Baseline-korrigierter Halb-Exponential-Fit fuer robustere Startwerte.
n_baseline = max(20, int(0.1 * len(Peak_shape)))
baseline = np.mean(Peak_shape[:n_baseline, 1])
y_corr = Peak_shape[:, 1] - baseline
A_guess = np.max(y_corr)
x0_guess = Peak_shape[y_corr > 0.5 * A_guess, 0].min()

tail_mask = (Peak_shape[:, 0] >= x0_guess) & (y_corr > 0.1 * A_guess)
tail_x = Peak_shape[tail_mask, 0]
tail_ln = np.log(y_corr[tail_mask])
l_guess = -np.polyfit(tail_x, tail_ln, 1)[0]

idx_fehler_grenze = np.argmin(np.abs(Peak_shape[:, 0] - 0.00123))
Fehler = np.std(Peak_shape[:idx_fehler_grenze, 1], ddof = 1)
p0 = [x0_guess, l_guess, A_guess]
popt, pcov = curve_fit(
    peak_fit,
    Peak_shape[:, 0],
    y_corr,
    p0=p0,
    bounds=([Peak_shape[:, 0].min(), 0.0, 0.0], [Peak_shape[:, 0].max(), np.inf, np.inf]),
    maxfev=50000,
    sigma = Fehler
)
fig, ax = plt.subplots(2, 1,figsize=(8,6), dpi=300, gridspec_kw={'height_ratios': (3,1), 'hspace': 0.05})

x_fit = np.linspace(Peak_shape[:, 0].min(), Peak_shape[:,  0].max(), 1000)
y_fit = peak_fit(x_fit, *popt) 
ax[0].plot(Peak_shape[:, 0], Peak_shape[:, 1]-baseline, 'o', label='Daten')
ax[0].plot(x_fit, y_fit, label='Peak-Fit', color='red')
ax[0].set_xlabel('Zeit [s]')
ax[0].set_ylabel('Amplitude [V]')
ax[0].legend()
ax[0].grid()
ax[0].text(0.6, 0.8, f'Fit-Parameter:\n$x_0$ = {popt[0]:.2e}$\pm$ {np.sqrt(pcov[0, 0]):.2e} s\n$\lambda$ = {popt[1]:.2e}$\pm$ {np.sqrt(pcov[1, 1]):.2e} 1/s\n$H$ = {popt[2]:.3f}$\pm$ {np.sqrt(pcov[2, 2]):.3f} V', transform=ax[0].transAxes, verticalalignment='top', fontsize=10,bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
residuals = Peak_shape[:, 1] - (peak_fit(Peak_shape[:, 0], *popt) + baseline)
ax[1].plot(Peak_shape[:,  0], residuals, 'o', color='maroon', label='Residuen', markersize=3)
ax[1].text(0.6, 0.4, r'$\frac{\chi^2}{N_{dof}}$ ' + f'= {np.sum((residuals / Fehler) ** 2)/(len(Peak_shape)-len(popt)):.2f}' + f', Ndof={len(Peak_shape)-len(popt)}', transform=ax[1].transAxes, verticalalignment='top', fontsize=12,bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))   
ax[1].axhline(0, color='black', linestyle='--', linewidth=1.5)
ax[1].set_xlabel('Zeit [s]')
ax[1].set_ylabel('Residuen')
ax[1].grid()
Flaeche = np.trapezoid(y_fit, x_fit)
plt.savefig('Peak_Fit.svg', dpi=300)
f = 1056
T = 20
Anzahl_perioden = f*T
P = 1#mikrowatt
P = ufloat(94, 2)
Energie_pro_Peak = P / f
A = ufloat(popt[2], np.sqrt(pcov[2, 2]))
print(f"Energie pro Peak: {Energie_pro_Peak:.3f} muW")
Spitzen_leistung = A*Energie_pro_Peak / Flaeche
print(f"Spitzenleistung: {Spitzen_leistung:.3f} muW")
print(Flaeche, popt[2])
# %%
