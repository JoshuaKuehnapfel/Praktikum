#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sympy import ln
from scipy.optimize import brentq
data_ = np.loadtxt('Helium_daten', dtype=str)
data = np.loadtxt('Daten', dtype=str)

data = np.append(data, data_, axis=0)

data = np.char.replace(data, ',', '.')

data = data.astype(float)

data_m = np.zeros((len(data)//3, data.shape[1]))
data_s = np.zeros((len(data)//3, data.shape[1]))

for i in range(len(data)//3):
    data_m[i,:] = data[3*i+1,:]
    data_s[i,:] = np.abs(data[3*i,:] - data[3*i+2,:])/2

data_C = data_m[:,1]

data_R = np.array([data_C[-1], 1361.74639, data_C[0]])

data_T = np.array([4.2, 77.5, 19.4+273.15])
ref_1 = data_m[0,0]
ref_2 = data_m[-1,0]
m = (77- (19.4+273.15))/(ref_2-ref_1)
b = 19.4+273.15 - m*ref_1
x_data = data_m[:,0]*m + b




plt.errorbar(x_data, data_m[:,2], yerr=data_s[:,2], fmt='.')
plt.xlabel('Temperatur [K]')
plt.ylabel('Widerstand [Ω]')
plt.title('Messwerte')
plt.show()  
print(x_data)

Silikon = data_m[:,4]
Silikon_s = data_s[:,4]
idx = np.where((data_m[:,4] <= 0.15*1e8))
Silikon = Silikon[idx]
Silikon_s = Silikon_s[idx]
x_data = x_data[idx]
plt.errorbar(x_data, Silikon, yerr=Silikon_s, fmt='.')
plt.grid()
plt.show()
plt.errorbar(x_data, Silikon, yerr=Silikon_s, fmt='.')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.show()
plt.errorbar(1/x_data,  Silikon, yerr=Silikon_s, fmt='.')
plt.yscale('log')
plt.grid()
plt.show()




def fit(x1, x2, x3, y1, y2, y3):
    assert(x1 < x2 < x3)
    assert(y1 > y2 > y3)

    d2 = x2 - x1
    d3 = x3 - x1
    L = (y2 - y1)/(y3 - y1)

    f = lambda c: (-np.expm1(-c*d2))/(-np.expm1(-c*d3)) - L
    c = brentq(f, 1e-12, 1.0)

    b = (y2 - y1)/(np.exp(-c*x2) - np.exp(-c*x1))
    a = y1 - b*np.exp(-c*x1)
    return a, b, c

a, b, c = fit(4.2, 77.5, 19.4+273.15, 2778.66, 1361.75, 1044.28)
x_fit = np.linspace(4.2, 19.4+273.15, 100)
y_fit = a + b*np.exp(-c*x_fit)
plt.plot(data_T, data_R, 'o', label='data points')
plt.plot(x_fit, y_fit, 'r-', label='Fit')
plt.legend()
plt.show()