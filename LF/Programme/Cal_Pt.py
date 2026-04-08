#%% 
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

Rausch = np.loadtxt('Rauschmessung', dtype=str)
Rausch = np.char.replace(Rausch, ',', '.')
Rausch = Rausch.astype(float)

sigma_Pt = np.std(Rausch[:,0] ,ddof=1)

data = np.loadtxt('Datei_Stickstoff', dtype=str)
data = np.char.replace(data, ',', '.')
data = data.astype(float)

data_m = np.zeros((len(data)//3, data.shape[1]))
data_s = np.zeros((len(data)//3, data.shape[1]))

for i in range(len(data)//3):
    data_m[i,:] = data[3*i+1,:]
    data_s[i,:] = np.std([data[3*i,:], data[3*i+1,:], data[3*i+2,:]], axis=0, ddof=1)
print(data_s[:,0], sigma_Pt)
data_m = data
raumtemp = 19.4 + 273.15
ref_1 = data_m[0,0]
ref_2 = data_m[-1,0]
m = (77.35 - raumtemp)/(ref_2-ref_1)
b = raumtemp - m*ref_1

sava_parameters = np.array([m, b])
np.savetxt('parameters_Pt', sava_parameters, delimiter=' ', fmt='%.5f')

x_data = data_m[:,0]*m + b

data_Pt = data_m[:,0] 
data_Pt_err = data_s[:,0] 

Fit_data_y = np.array([data_Pt[0], data_Pt[-1]])
Fit_data_x = np.array([292.55, 77.35])

y = np.linspace(15, 115, 39)

lit = np.linspace(50, 300, 100)

plt.errorbar(Fit_data_x, Fit_data_y, yerr = sigma_Pt, fmt='o', label='data points (fit)')
plt.plot(m*y+b, y, label='linear fit')
plt.errorbar(m*data_Pt +b, data_Pt, yerr=sigma_Pt, marker='.', linestyle='', label='measured data points transformed', zorder=3)
plt.plot(lit, 39.74/100 * lit -10.55, label='literature values', zorder=2, linestyle='dashed')
plt.ylabel(r'$R$ / $\Omega$')
plt.xlabel(r'$T$ / K')
plt.title('Linear fit of the Pt resistance as a function of the temperature')
plt.legend()
plt.grid()
plt.show()

print("bei T=0°C Widerstand: Fit  =", (273.15-b)/m, "\n Literatur = 100", "Differenz = ", (273.15-b)/m - 100)

print("bei T=-100°C Widerstand: Fit  =", (173.15-b)/m, "\n Literatur = 60.26", "Differenz = ", (173.15-b)/m - 60.26)

#Quelle: https://www.sab-kabel.de/kabel-konfektion-temperaturmesstechnik/technische-daten/temperaturmesstechnik/mantel-widerstandsthermometer/pt100-tabelle-widerstandstabelle-messwiderstaende.html?srsltid=AfmBOopi4tQzif3Bn6-c2lWsd-SDS6QpQtRW-yIFq5gOKeP06DnZYh4i

T = m*data_Pt + b

print(T)

save_data = np.column_stack((T))
np.savetxt('T_daten_hoch', save_data, delimiter=' ', fmt='%.10f')

# %%
