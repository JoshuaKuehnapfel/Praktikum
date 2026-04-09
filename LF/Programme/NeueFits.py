#%%
from Funktionen import *
from scipy.optimize import brentq
data_R = np.array([data_H_m[-1, 1], data_N_m[-1, 1], data_N_m[0, 1]])
data_T = np.array([4.22, 77.35, 18.4+273.15])

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

popt, pcov = plot_fit(f, data_T, data_R, np.ones_like(data_T)*0.0001, sigma_C, p0)
ax = plt.gca()

plt.title('Exponential fit of the carbon resistance as a function of temperature')
plt.text(0.6, 0.67, r'Fit-function= $a + b \cdot e^{\frac{c}{T}}$', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
plt.text(0.6, 0.6, f'a = {popt[0]:.2f}$\pm$ {np.sqrt(pcov[0,0]):.2f}\nb = {popt[1]:.2f}$\pm$ {np.sqrt(pcov[1,1]):.2f}\nc = {popt[2]:.2f}$\pm$ {np.sqrt(pcov[2,2]):.2f}', transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')
ax.set_xlabel(r'$T$ / K')
ax.set_ylabel(r'$R$ / $\Omega$')


plt.savefig('Cal_C_fit.svg', dpi = 300)

#




