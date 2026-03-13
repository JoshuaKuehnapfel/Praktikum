
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.optimize import curve_fit
import scipy.constants as const
import os



plt.rcParams['font.size'] = 16.0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelsize'] = 'medium'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['lines.linewidth'] = 2.0

Fehlerkorrektur = 1
N = ufloat(322, 2.8)
e = const.e

U_R = np.array([3100, 3400, 3300, 3200, 3350, 3000, 3250, 3375, 3325])
U0 = np.array([2900, 3000, 3100, 3200, 3300, 3350, 3400])
U_80_20 = np.array([2800, 2880, 2960, 3040, 3120, 3210, 3300])
U_60_40 = np.array([2900, 3000, 3100, 3200, 3300, 3400, 3450, 3500])
U_G1T = np.array([2576, 2579, 2582, 2585, 2588, 2592, 2596])
R_an = np.array([2207, 14639, 13502, 7456, 14360, 2, 12109, 14294, 13801])
R_aus = np.array([3, 39, 20, 5, 30, 0, 1, 23, 19])


R = R_an - R_aus
R_uns = np.sqrt(R_an + R_aus)/60
sort_indices = np.argsort(U_R)
R = R[sort_indices]/60
R_uns = R_uns[sort_indices]
U_R = U_R[sort_indices]
Rate = np.average(R[6:], weights = 1/R_uns[6:]**2,axis=0)
Rate_uns = np.sqrt(1/np.sum(1/R_uns[6:]**2))

Daten_namen = ['V_Drift', 'V_Drift_80_20', 'V_Drift_60_40']


fig, ax = plt.subplots(1, 1, figsize=(10, 8))
ax.hlines(Rate,min(U_R),max(U_R), color='red', linestyle='dashed', label=f'Rate: {Rate:.2f} ± {Rate_uns:.2f} Hz')
plt.fill_between(U_R, Rate + Rate_uns, Rate - Rate_uns, alpha=0.3, color='lightcoral', label='1σ Unsicherheitsbereich')
ax.errorbar(U_R, R, yerr=R_uns*2, fmt='none', color = 'black', linewidth = 2, label = 'Unsicherheiten (2σ)')
ax.plot(U_R, R, linestyle='dashdot', color = 'black')
ax.scatter(U_R[:6], R[:6], c='red', label='Datenpunkte', zorder=5)
ax.scatter(U_R[6:], R[6:], c='darkblue', label='Datenpunkte Ratenbestimmung', zorder=5)

ax.set_xlabel(r'Spannung $V_{Drift}$ [V]')
ax.set_ylabel('Rate [Hz]')
ax.legend()
ax.grid()
ax.set_ylabel('Rate [Hz]')
ax.legend(loc = 'lower right')
plt.grid(True)
plt.tight_layout()
plt.savefig('Bilder/Zählraten.svg', dpi=300)
plt.clf()
plt.close()


def exp_fit (x, a, b):
    return a * np.exp(b * x)
p0 = {'V_Drift': (1.4e-7, 7e-3),
      'V_Drift_80_20': (1e-6, 1e-3),
      'V_Drift_60_40': (1e-6, 1e-3),
      'V_G1T': (1e-7, 1e-3)}
Legende_daten = {'V_Drift': r'$V_{Drift}$ bei 70/30 (Ar/CO2)',
                  'V_Drift_80_20': r'$V_{Drift}$ bei 80/20 (Ar/CO2)',
                  'V_Drift_60_40': r'$V_{Drift}$ bei 60/40 (Ar/CO2)',
                  'V_G1T': r'$V_{G1T}$'}
def plot_errorfit(Data, U, Gain, Gain_uns, c='orange'):
    if Data == 'V_Drift_60_40':
        Gain = Gain[3:]
        Gain_uns = Gain_uns[3:]
        U = U[3:]
        popt, pcov = curve_fit(exp_fit, U, Gain, sigma=Gain_uns, absolute_sigma=True, p0 = p0[Data])
        x = np.linspace(min(U), max(U)+20, 300)
        sa2 = pcov[0,0]
        sb2 = pcov[1,1]
        sab = pcov[0,1]
        a = popt[0]
        b = popt[1]
        y = exp_fit(x, *popt)
        sigma = np.sqrt(
            np.exp(2*b*x)*sa2
            + (a*x*np.exp(b*x))**2 * sb2
            + 2*np.exp(b*x)*(a*x*np.exp(b*x))*sab
        )
        plt.plot(x, y)
        plt.fill_between(x, y-sigma, y+sigma, alpha=0.3, color=c, label='1σ Unsicherheitsbereich')
        plt.plot(x, y, label='Exponentielle Anpassung', color=c)
        plt.title('Gain als Funktion von '+Legende_daten[Data])
    else:
        popt, pcov = curve_fit(exp_fit, U, Gain, sigma=Gain_uns, absolute_sigma=True, p0 = p0[Data])
        x = np.linspace(min(U), max(U)+20, 300)
        sa2 = pcov[0,0]
        sb2 = pcov[1,1]
        sab = pcov[0,1]
        a = popt[0]
        b = popt[1]
        y = exp_fit(x, *popt)
        sigma = np.sqrt(
            np.exp(2*b*x)*sa2
            + (a*x*np.exp(b*x))**2 * sb2
            + 2*np.exp(b*x)*(a*x*np.exp(b*x))*sab
        )
        plt.plot(x, y)
        plt.fill_between(x, y-sigma, y+sigma, alpha=0.3, color=c, label='1σ Unsicherheitsbereich')
        plt.plot(x, y, label='Exponentielle Anpassung', color=c)
        plt.title('Gain als Funktion von '+Legende_daten[Data])
Rate = ufloat(Rate, Rate_uns)

def plotting(data, U, c = 'darkblue'):
    Strom, Strom_err = np.loadtxt(f'Daten/{data}.csv', delimiter=',',  unpack=True, usecols=(1,2))
    # if data == 'V_Drift_60_40':
    #     Strom, Strom_err = np.loadtxt(f'Daten/{data}.csv', delimiter=',',  unpack=True, usecols=(1,2), skiprows = 6)
    #     U = U[3:]
    I = Strom[::2]
    I_err = Strom_err[::2]/Fehlerkorrektur
    I_aus = Strom[1::2]
    I_aus_err = Strom_err[1::2]/Fehlerkorrektur
    Gain = np.zeros(len(I))
    Gain_uns = np.zeros(len(I))

    for i in range(len(I)):
        I_temp = ufloat(I[i], I_err[i])
        I_temp_aus = ufloat(I_aus[i], I_aus_err[i])
        G = np.abs((I_temp-I_temp_aus)/(e*Rate*N))
        Gain[i] = G.n
        Gain_uns[i] = G.s
    plt.errorbar(U, Gain, yerr=Gain_uns, fmt='o', label='Gain mit Unsicherheiten', linestyle='dashdot', color=c)
    plt.ylabel('Gain')
    plt.grid()
    return Gain, Gain_uns, I, I_err, I_aus, I_aus_err

for data, U in zip(Daten_namen, [U0, U_80_20, U_60_40]):
    fig, ax = plt.subplots(3, 1, figsize=(10, 12), sharex=True, gridspec_kw={'height_ratios': [2, 1, 1]})
    plt.sca(ax[0])
    Gain, Gain_uns, I, I_err, I_aus, I_aus_err = plotting(data, U)
    plot_errorfit(data, U, Gain, Gain_uns)
    ax[1].scatter(U, -I_err/I*100, label='Rel. Fehler Strom (Quelle an)', color = 'maroon')
    ax[1].set_label('Strom I [A]')
    ax[1].set_ylabel('Rel. Fehler I [%]')
    ax[2].set_xlabel(r'Spannung $V_{Drift}$ [V]')
    ax[2].scatter(U, -I_aus_err/I_aus*100, label='Rel. Fehler Strom (Quelle aus)', color = 'maroon')
    ax[2].set_ylabel(r'Rel. Fehler $I_{aus}$ [%]')
    for a in ax:
        a.legend(loc = 'upper left')
        a.grid(True)
    plt.savefig(f'Bilder/Verstärkung_{data}.svg', dpi=300)


