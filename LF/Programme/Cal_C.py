#%% 
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.optimize import brentq

data_H = np.loadtxt('Helium_daten', dtype=str)
data_N = np.loadtxt('Datei_Stickstoff', dtype=str)

data_H = np.char.replace(data_H, ',', '.')
data_N = np.char.replace(data_N, ',', '.')

data_H = data_H.astype(float)
data_N = data_N.astype(float)

data_m = np.zeros((len(data_H)//3, data_H.shape[1]))
data_s = np.zeros((len(data_H)//3, data_H.shape[1]))

for i in range(len(data_H)//3):
    data_m[i,:] = data_H[3*i+1,:]
    data_s[i,:] = np.abs(data_H[3*i,:] - data_H[3*i+2,:])/2
data_m = data_H
data_C_H = data_m[:,1]

data_m_ = np.zeros((len(data_N)//3, data_N.shape[1]))
data_s_ = np.zeros((len(data_N)//3, data_N.shape[1]))

for i in range(len(data_N)//3):
    data_m_[i,:] = data_N[3*i+1,:]
    data_s_[i,:] = np.abs(data_N[3*i,:] - data_N[3*i+2,:])/2

data_C_N = data_m_[:,1]

data_R = np.array([data_C_H[-1], data_C_N[-1], data_C_N[0]])

data_T = np.array([4.22, 77.35, 18.4+273.15])

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

a, b, c = fit(*data_T, *data_R)
x_fit = np.linspace(0.1, 20+273.15, 300)
y_fit = a + b*np.exp(c/x_fit)
plt.plot(data_T, data_R, 'o', label='data points')
plt.plot(x_fit, y_fit, 'r-', label='Fit')
plt.title('Exponential fit of the carbon resistance as a function of temperature')
plt.xlabel(r'$T$ / K')
plt.ylabel(r'$R$ / Ω')
plt.grid()
plt.legend()
plt.show()
T_niedrig = c / (np.log((data_C_H - a)/b))
print(a, b, c)

#print("Vergleich Literatur Kohlewiderstand: T bei R=100Ω = ", (-1/c)*np.log((100 - a)/b), "K", "\n Literatur = ")

save_data = np.column_stack((T_niedrig))
np.savetxt('T_daten_niedrig', save_data, delimiter=' ', fmt='%.2f')
print(T_niedrig)
#Quelle Stickstoff: https://www.chemie.de/lexikon/Stickstoff.html
#Quelle Helium: https://www.chemie.de/lexikon/Helium.html 