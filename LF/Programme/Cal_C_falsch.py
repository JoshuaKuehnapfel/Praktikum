#%% 
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.optimize import brentq

data_ = np.loadtxt('Helium_daten', dtype=str)
data = np.loadtxt('Datei_Stickstoff', dtype=str)

data = np.append(data, data_, axis=0)


data_ = np.char.replace(data_, ',', '.')

data_ = data_.astype(float)

data = np.char.replace(data, ',', '.')

data = data.astype(float)

data_m = np.zeros((len(data)//3, data.shape[1]))
data_s = np.zeros((len(data)//3, data.shape[1]))

for i in range(len(data)//3):
    data_m[i,:] = data[3*i+1,:]
    data_s[i,:] = np.abs(data[3*i,:] - data[3*i+2,:])/2

data_C = data_m[:,1]

data_m_ = np.zeros((len(data_)//3, data_.shape[1]))
data_s_ = np.zeros((len(data_)//3, data_.shape[1]))

for i in range(len(data_)//3):
    data_m_[i,:] = data_[3*i+1,:]
    data_s_[i,:] = np.abs(data_[3*i,:] - data_[3*i+2,:])/2

data_C_ = data_m_[:,1]


'''
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

# Jetzt hast du a, b und dein gewähltes c.
# Diese Werte nutzt du nun als STARTWERTE für den echten exp_fun Fit.
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

plt.legend()
plt.grid()
plt.show()
'''

data_R = np.array([data_C[-1], 1361.74639, data_C[0]])

data_T = np.array([4.2, 77.5, 19.4+273.15])

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

a, b, c = fit(data_T[0], data_T[1], data_T[2], data_R[0], data_R[1], data_R[2])
x_fit = np.linspace(4.2, 19.4+273.15, 100)
y_fit = a + b*np.exp(-c*x_fit)
plt.plot(data_T, data_R, 'o', label='data points')
plt.plot(x_fit, y_fit, 'r-', label='Fit')
plt.title('Exponential fit of the carbon resistance as a function of temperature')
plt.xlabel(r'$T$ / K')
plt.ylabel(r'$R$ / Ω')
plt.grid()
plt.legend()
plt.show()

T = (-1/c)*np.log((data_C_ - a)/b)

print("Vergleich Literatur Kohlewiderstand: T bei R=100Ω = ", (-1/c)*np.log((100 - a)/b), "K", "\n Literatur = ")



save_data = np.column_stack((T))
np.savetxt('T_daten_niedrig', save_data, delimiter=' ', fmt='%.2f')