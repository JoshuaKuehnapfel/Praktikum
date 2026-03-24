
#%%
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
from scipy.optimize import curve_fit
daten = np.array([2.97, 2.94, 2.97, 2.87])
daten_uns = np.array([0.05, 0.05, 0.05, 0.1])

weights = 1/daten_uns**2
mean = np.average(daten, weights=weights)
std_dev = np.sqrt(1/np.sum(weights))
# print(f"Gewichteter Mittelwert: {mean:.2f} +- {std_dev:.2f}")

rho = 1.2041*1e-3

d = ufloat(2.95, 0.12)

Energie = 4.69
stopping_power = ufloat(700, 10)
delta_Energie = Energie-stopping_power*rho * d
print(f"E_0 = {Energie} MeV, E_detektiert: {delta_Energie:.2f} MeV")
print('-' * 30)

Energie = 5.3
stopping_power = ufloat(750, 10)
delta_Energie = Energie-stopping_power*rho * d
print(f"E_0 = {Energie} MeV, E_detektiert: {delta_Energie:.2f} MeV")
print('-' * 30)

Energie = 5.69
stopping_power = ufloat(760, 10)
delta_Energie = Energie-stopping_power*rho * d
print(f"E_0 = {Energie} MeV, E_detektiert: {delta_Energie:.2f} MeV")
print('-' * 30)

Energie = 6.00
stopping_power = ufloat(799.8, 0.1)
delta_Energie = Energie-stopping_power*rho * d
print(f"E_0 = {Energie} MeV, E_detektiert: {delta_Energie:.2f} MeV")
print('-' * 30)
delta_Energien = np.zeros(4)
delta_energien_unc = np.zeros(4)
abstände = np.array([0, 0.99, 1.79, 4.19])
abstände_uns = np.array([0.05, 0.05, 0.05, 0.1])
Energie = 7.69
stopping_power = ufloat(899.8, 0.1)
delta_Energie = Energie-stopping_power*rho * d
delta_Energien[0] = delta_Energie.n
delta_energien_unc[0] = delta_Energie.s
print(f"E_0 = {Energie} MeV, E_detektiert: {delta_Energie:.2f} MeV")
print('-' * 30)

print("Energieverlust für 7.69 MeV mit d = 2.95 cm:")

d_n =  ufloat(0.99, 0.05)
delta_Energie = Energie-stopping_power*rho * (d_n+d)
delta_Energien[1] = delta_Energie.n
delta_energien_unc[1] = delta_Energie.s
print(f"E_0 = {Energie} MeV, E_detektiert({d_n}): {delta_Energie:.2f} MeV")
print('-' * 30)
d_n =  ufloat(1.79, 0.05)
delta_Energie = Energie-stopping_power*rho * (d_n+d)
delta_Energien[2] = delta_Energie.n
delta_energien_unc[2] = delta_Energie.s
print(f"E_0 = {Energie} MeV, E_detektiert({d_n}): {delta_Energie:.2f} MeV")
print('-' * 30)
d_n =  ufloat(4.19, 0.1)
delta_Energie = Energie-stopping_power*rho * (d_n+d)
delta_Energien[3] = delta_Energie.n
delta_energien_unc[3] = delta_Energie.s
print(f"E_0 = {Energie} MeV, E_detektiert({d_n}): {delta_Energie:.2f} MeV")
print('-' * 30)

def linear(x, m, b):
    return m * x + b

abstände = np.array([0, 0.99, 1.79, 4.19])
abstände_uns = np.array([0.05, 0.05, 0.05, 0.1])
delta_Energien = np.array([7.69, 7.5, 6.5, 5.5])
delta_Energien_unc = np.array([0.5, 0.5, 0.5, 0.5])

popt, pcov = curve_fit(linear, delta_Energien, abstände, sigma=np.sqrt(np.array(delta_Energien_unc)**2 + np.array(abstände_uns)**2), absolute_sigma=True)
x_smooth = np.linspace(min(delta_Energien), max(delta_Energien), 100)
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8,6), dpi=300, gridspec_kw={'height_ratios': (3,1), 'hspace': 0.05})
ax[0].errorbar(delta_Energien, abstände, xerr = delta_Energien_unc, yerr = abstände_uns, fmt = 'o', label='Energie $E_{detektiert}$', zorder=2)
ax[0].plot(x_smooth, linear(x_smooth, *popt), color='black', label='Fit', zorder=2)
ax[0].set_xlabel(r'Energie $E_{detektiert}$ (MeV)')
ax[0].set_ylabel('Entfernung vom Referenzpunkt (cm)')
ax[0].text(0.05, 0.45, f'm = {popt[0]:.2f} +- {np.sqrt(pcov[0,0]):.2f}\nb = {popt[1]:.2f} +- {np.sqrt(pcov[1,1]):.2f}', transform=ax[0].transAxes, verticalalignment='top', fontsize = 12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
ax[1].errorbar(delta_Energien, abstände - linear(delta_Energien, *popt), xerr=delta_Energien_unc, yerr=abstände_uns, fmt='o', color='maroon', label='Residuen')
ax[1].axhline(0, linewidth=1.5, linestyle = 'dotted',color='black')
ax[1].set_xlabel(r'Energie $E_{detektiert}$ (MeV)')
ax[1].set_ylabel('Residuen')
ax[1].text(0.25, 0.95, r'Reduziertes $\frac{\chi^2}{N_{dof}}$:'+ f' {np.sum((abstände - linear(delta_Energien, *popt))**2 / (abstände_uns**2)) / (len(delta_Energien) - 2):.3f}' + f'; Ndof = {len(delta_Energien) - 2}', transform=ax[1].transAxes, verticalalignment='top', fontsize=12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
for a in ax:
    a.grid()
    a.legend()
fig.suptitle('Kalibrierung mit linearem Fit')
plt.savefig('Kalibrierung_mit_Fit_7.69.svg', dpi=300)


m = ufloat(-790.2, 83.6)
b = ufloat(76581,44.5)
M = ufloat(-1.72, 0.29)
B = ufloat(13.43, 1.99)

m_prime = m * M
b_prime = m * B + b
print(f"m' = {m_prime:.2f} +- {m_prime.s:.2f}, b' = {b_prime:.2f} +- {b_prime.s:.2f}")
print('-'*30)
abstände = np.array([0, 0.99, 1.79, 4.19])
abstände_uns = np.array([0.00, 0.05, 0.05, 0.1])
RLit = ufloat(7.06, 0.35)
for i in range(4):
    R = ufloat(abstände[i], abstände_uns[i])
    R_rest = RLit - R
    print(f"Abstand: {R:.2f} cm, Restreichweite: {R_rest:.2f}")