#%% 
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

data = np.loadtxt('Datei_Stickstoff', dtype=str)

data = np.char.replace(data, ',', '.')

data = data.astype(float)

data_m = np.zeros((len(data)//3, data.shape[1]))
data_s = np.zeros((len(data)//3, data.shape[1]))

for i in range(len(data)//3):
    data_m[i,:] = data[3*i+1,:]
    data_s[i,:] = np.abs(data[3*i,:] - data[3*i+2,:])/2
ref_1 = data_m[0,0]
ref_2 = data_m[-1,0]
m = (77- (19.4+273.15))/(ref_2-ref_1)
b = 19.4+273.15 - m*ref_1
x_data = data_m[:,0]*m + b

data_Pt = data_m[:,0] 
data_Pt_err = data_s[:,0] 

data_C = data_m[:,1]
data_C_err = data_s[:,1]

data_Cu = data_m[:,2] 
data_Cu_err = data_s[:,2]

data_Ta = data_m[:,3]
data_Ta_err = data_s[:,3]

data_Sz = data_m[:32,4]
data_Sz_err = data_s[:32,4]

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
    return popt, pcov

plt.plot(x_data, data_Pt, linestyle='--', marker='.')
plt.xlabel(r'$T$ / K')
plt.ylabel(r'$R / Ω $')
plt.grid()


popt, pcov= fit_mit_plot(linear, x_data, data_Pt, data_Pt_err)
print(popt[0])



print(data_Pt[0], max(data_Pt_err))
print(data_Pt[-1], data_Pt_err[-1])

print(m, b)
 
print(data_C[0], data_C[-1])
# %%
