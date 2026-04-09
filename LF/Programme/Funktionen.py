#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData
from scipy.optimize import curve_fit
import os
from scipy.optimize import brentq
import uncertainties as unc
import uncertainties.unumpy as unp
os.chdir('G:/Uni/Praktikum/LF/Programme')


Rausch = np.loadtxt('Rauschmessung', dtype=str)
Rausch = np.char.replace(Rausch, ',', '.')
Rausch = Rausch.astype(float)
sigma_Pt = np.std(Rausch[:,0] ,ddof=1)
sigma_C = np.std(Rausch[:,1] ,ddof=1)
data_N = np.loadtxt('Datei_Stickstoff', dtype=str)
data_N = np.char.replace(data_N, ',', '.')
data_N = data_N.astype(float)

data_N_m = np.zeros((len(data_N)//3, data_N.shape[1]))
data_N_s = np.zeros((len(data_N)//3, data_N.shape[1]))

data_H = np.loadtxt('Helium_daten', dtype=str)
data_H = np.char.replace(data_H, ',', '.')
data_H = data_H.astype(float)

data_H_m = np.zeros((len(data_H)//3, data_H.shape[1]))
data_H_s = np.zeros((len(data_H)//3, data_H.shape[1]))

for i in range(len(data_N)//3):
    data_N_m[i,:] = data_N[3*i+1,:]
    data_N_s[i,:] = 1/3*np.abs([data_N[3*i,:]+ data_N[3*i+2,:]-2*data_N[3*i+1,:]])
for i in range(len(data_H)//3):   
    data_H_m[i,:] = data_H[3*i+1,:]
    data_H_s[i,:] = 1/3*np.abs([data_H[3*i,:]+ data_H[3*i+2,:]-2*data_H[3*i+1,:]])


data_R = np.array([data_H_m[-1, 1], data_N_m[-1, 1], data_N_m[0, 1]])
data_T = np.array([4.22, 77.35, 18.4+273.15])


def linear_func(p, x):
    m, c = p
    return m * x + c

def fit_xerr_yerr(func, xdata, ydata, xerr, yerr, p0):
    data = RealData(xdata, ydata, sx=xerr, sy=yerr)
    model = Model(func)
    odr = ODR(data, model, beta0=p0)
    output = odr.run()
    return output.beta, output.cov_beta

def plot_fit(func, xdata, ydata, xerr, yerr, p0):
    popt, pcov = fit_xerr_yerr(func, xdata, ydata, xerr, yerr, p0)
    a_fit, b_fit, c_fit = popt
    print(f"Fitted parameters: a = {a_fit:.5f}, b = {b_fit:.5f}, c = {c_fit:.5f}")
    print(f"Covariance matrix:\n{pcov}")
    # Plotting the data and the fit
    plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, fmt='o', label='Data with error bars')
    x_fit = np.linspace(min(xdata), max(xdata), 100)
    y_fit = func(popt, x_fit)
    plt.plot(x_fit, y_fit, label='ODR Fit', color='red')
    plt.legend()
    plt.grid()
    return popt, pcov



#Kohlenstoff Exponential fit

def fit(x1, x2, x3, y1, y2, y3):
    # Annahme: x sind Temperaturen (T), y sind Widerstände (R)
    assert(x1 < x2 < x3)
    assert(y1 > y2 > y3)

    # Hier liegt die Änderung: Differenzen der Kehrwerte für das 1/T Modell
    inv_d2 = 1/x2 - 1/x1

    inv_d3 = 1/x3 - 1/x1

    L = (y2 - y1)/(y3 - y1)

    f = lambda c: np.expm1(c * inv_d2) / np.expm1(c * inv_d3) - L
    c = brentq(f, -1000, 5000)

    b = (y2 - y1) / (np.exp(c/x2) - np.exp(c/x1))
    a = y1 - b * np.exp(c/x1)
    return a, b, c
p0 = fit(*data_T, *data_R)

def f(p, x):
    a, b, c = p
    return a + b * np.exp(c / x)

popt, pcov = plot_fit(f, data_T, data_R, np.ones_like(data_T)*1e-30, sigma_C, p0)
ax = plt.gca()

plt.title('Exponential fit of the carbon resistance as a function of temperature')
plt.text(0.6, 0.67, r'Fit-function= $a + b \cdot e^{\frac{c}{T}}$', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
plt.text(0.6, 0.6, f'a = {popt[0]:.2f}$\pm$ {np.sqrt(pcov[0,0]):.2f}\nb = {popt[1]:.2f}$\pm$ {np.sqrt(pcov[1,1]):.2f}\nc = {popt[2]:.2f}$\pm$ {np.sqrt(pcov[2,2]):.2f}', transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')
ax.set_xlabel(r'$T$ / K')
ax.set_ylabel(r'$R$ / $\Omega$')


plt.savefig('Cal_C_fit.svg', dpi = 300)

data_C_H = unp.uarray(data_H_m[:,1], data_H_s[:,1])

a, b, c = unc.correlated_values(popt, pcov)
Temperatur_niedrig = c / (unp.log((data_C_H - a)/b))
Temperatur = np.column_stack((unp.nominal_values(Temperatur_niedrig), unp.std_devs(Temperatur_niedrig)))
print(Temperatur)
np.savetxt('T_daten_niedrig_mit_uns', Temperatur, delimiter=' ', fmt='%.4e')

raumtemp = 19.4 + 273.15
ref_1 = unc.ufloat(data_N_m[0,0], sigma_Pt)
ref_2 = unc.ufloat(data_N_m[-1,0], sigma_Pt)
m = (77.35 - raumtemp)/(ref_2-ref_1)
b = raumtemp - m*ref_1



data_pt_N = unp.uarray(data_N_m[:,0], data_N_s[:,0])
Temperatur_hoch = m*data_pt_N + b
data_Pt = data_N_m[:,0]
save_data = np.column_stack((unp.nominal_values(Temperatur_hoch), unp.std_devs(Temperatur_hoch)))
np.savetxt('T_daten_hoch_mit_uns', save_data, delimiter=' ', fmt='%.4e')

Fit_data_y = np.array([data_Pt[0], data_Pt[-1]])
Fit_data_x = np.array([292.55, 77.35])

plt.clf()
y = np.linspace(15, 115, 39)

lit = np.linspace(50, 300, 100)

plt.errorbar(Fit_data_x, Fit_data_y, yerr = sigma_Pt, fmt='o', label='data points (fit)')
plt.plot(m.n*y+b.n, y, label='linear fit')
plt.errorbar(unp.nominal_values(Temperatur_hoch), data_Pt, xerr = unp.std_devs(Temperatur_hoch), yerr=sigma_Pt, marker='.', linestyle='', label='measured data points transformed', zorder=3)
plt.plot(lit, 39.74/100 * lit -10.55, label='literature values', zorder=2, linestyle='dashed')
plt.ylabel(r'$R$ / $\Omega$')
plt.xlabel(r'$T$ / K')
plt.text(0.5, 0.47, r'Fit-function: R(T) = $m \cdot T + b$', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
plt.text(0.55, 0.4, f'm = {m.n:.5f}$\pm$ {m.s:.5f}\nb = {b.n:.3f}$\pm$ {b.s:.3f}', transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')
plt.title('Linear fit of the Pt resistance as a function of the temperature')
plt.legend()
plt.grid()
plt.show()

plt.clf()
plt.figure(figsize=(8,5), layout='tight')
plt.plot(unp.std_devs(Temperatur_niedrig)/unp.nominal_values(Temperatur_niedrig)*100, marker='o', linestyle='', label='Carbon-scale', zorder=3)
plt.plot(unp.std_devs(Temperatur_hoch)/unp.nominal_values(Temperatur_hoch)*100, marker='o', linestyle='', label='Platinum-scale', zorder=3)
plt.grid()
plt.xlabel('data point index')
plt.ylabel('relative uncertainty in [%]')
plt.legend()
print(sigma_C/data_H_m[-1,1], sigma_Pt/data_H_m[-1,0])
plt.title('Relative uncertainties of the temperature values \n obtained from the Pt and C resistance measurements')
plt.savefig('relative_uncertainties.svg', dpi=300)