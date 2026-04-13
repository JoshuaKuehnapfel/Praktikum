import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

data_ = np.loadtxt('Helium_daten', dtype=str)
data = np.loadtxt('Datei_Stickstoff', dtype=str)

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

def exp_fun(T, a, b, c):
    return a * np.exp(b / T) + c

def linear(T, m, n):
    return m * T + n

c_guess = np.min(data_R) * 0.95

y_linear = np.log(data_R - c_guess)
x_linear = 1 / data_T

popt_l, cov_l = curve_fit(linear, x_linear, y_linear)
b_fit = popt_l[0]
a_fit = np.exp(popt_l[1])
p0 = [a_fit, b_fit, c_guess]

popt_l, cov_l = curve_fit(linear, np.log(data_T), np.log(data_R))

def exp_fun(T, a, b, c):
    return a * np.exp(b / T) + c

p0 = [300, 20, 900]

popt, pcov = curve_fit(exp_fun, data_T, data_R, p0=p0, bounds=(0, [np.inf, np.inf, np.inf]))

x_ = np.linspace(3, 300, 1000)
y_ = exp_fun(x_, *popt)

plt.plot(data_T, data_R, 'o', label='data points')
#plt.plot(popt[1]/ (np.log((data_C - popt[2]) / popt[0] )), data_C, marker='.', linestyle='', label='data points transformed', zorder=3)
plt.plot(x_, y_, label='exponential fit')
plt.xlabel(r'$T$ / K')
plt.ylabel(r'$R$ / Ω')
plt.title('Exponential fit of the carbon resistance as a function of temperature')
plt.legend()
plt.grid()
plt.show()