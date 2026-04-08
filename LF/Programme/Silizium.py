#%% 
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

T_data_N = np.loadtxt('T_daten_hoch', dtype=str)
T_data_H = np.loadtxt('T_daten_niedrig', dtype=str)

T_data_N = T_data_N.astype(float)
T_data_H = T_data_H.astype(float)

T_data_N = T_data_N[:33]
T_data_H = T_data_H[:16]

#T_data = np.append(T_data_1, T_data_2)
T_data_N_log = np.log(T_data_N.astype(float))
T_data_H_log = np.log(T_data_H.astype(float))

data_N = np.loadtxt('Datei_Stickstoff', dtype=str)
data_H = np.loadtxt('Helium_daten', dtype=str)

data_N = np.char.replace(data_N, ',', '.')
data_H = np.char.replace(data_H, ',', '.')

data_N = data_N.astype(float)
data_H = data_H.astype(float)


data_m = np.zeros((len(data_N)//3, data_N.shape[1]))
data_s = np.zeros((len(data_N)//3, data_N.shape[1]))

for i in range(len(data_N)//3):
    data_m[i,:] = data_N[3*i+1,:]
    data_s[i,:] = np.abs(data_N[3*i,:] - data_N[3*i+2,:])/2

data_Sz_N = data_m[:,4] 
data_Sz_N = data_Sz_N[:33] #nutzbare Werte
#data_Sz_err_N = data_s[:,4]

data_m_ = np.zeros((len(data_H)//3, data_H.shape[1]))
data_s_ = np.zeros((len(data_H)//3, data_H.shape[1]))


for i in range(len(data_H)//3):
    data_m_[i,:] = data_H[3*i+1,:]
    data_s_[i,:] = np.abs(data_H[3*i,:] - data_H[3*i+2,:])/2
    
data_Sz_H = data_m_[:,4] 
data_Sz_H = data_Sz_H[:16] #nutzbare Werte
#data_Sz_err_H = data_s_[:,4]



plt.plot(T_data_N, data_Sz_N, marker='.', linestyle='', label='data points (pt temp. scale)', zorder=3)
plt.plot(T_data_H, data_Sz_H, marker='.', linestyle='', label='data points (c temp. scale)', zorder=3)
plt.xlabel(r'$T$ / K')
plt.ylabel(r'$R$ / $\Omega$')
plt.title('Plot of the silicon resistance as a function of the temperature')
plt.legend()
plt.grid()
plt.show()


plt.plot(-T_data_N_log, -np.log(data_Sz_N), marker='.', linestyle='', label='data points (pt temp. scale)', zorder=3)
plt.plot(-T_data_H_log, -np.log(data_Sz_H), marker='.', linestyle='', label='data points (c temp. scale)', zorder=3)
plt.xlabel(r'$ln(1/T)$')
plt.ylabel(r'$ln(1/R)$')
plt.title('Plot of the silicon resistance as a function of the temperature')
plt.legend()
plt.grid()
plt.show()


plt.plot(T_data_N_log, -np.log(data_Sz_N), marker='.', linestyle='', label='data points (pt temp. scale)', zorder=3)
plt.plot(T_data_H_log, -np.log(data_Sz_H), marker='.', linestyle='', label='data points (c temp. scale)', zorder=3)
plt.xlabel(r'$ln(T)$')
plt.ylabel(r'$ln(1/R)$')
plt.title('Plot of the silicon resistance as a function of the temperature')
plt.legend()
plt.grid()
plt.show()

plt.plot(1/T_data_N_log, np.log(data_Sz_N), marker='.', linestyle='', label='data points (pt temp. scale)', zorder=3)
plt.plot(1/T_data_H_log, np.log(data_Sz_H), marker='.', linestyle='', label='data points (c temp. scale)', zorder=3)
plt.xlabel(r'$1/T$ / 1/K')
plt.ylabel(r'$ln(R)$')
plt.title('Plot of the silicon resistance as a function of the temperature')
plt.legend()
plt.grid()
plt.show()