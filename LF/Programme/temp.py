#%% 
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
import scipy.odr as odr
import os
#os.chdir('G:/Uni/Praktikum/LF/Programme')
plt.rcParams['font.size'] = 24.0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelsize'] = 'medium'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = 15

Rausch = np.loadtxt('Rauschmessung', dtype=str)
Rausch = np.char.replace(Rausch, ',', '.')
Rausch = Rausch.astype(float)

sigma_Cu = np.abs(np.std(Rausch[:,2] ,ddof=1))


par = np.loadtxt('parameters_Pt', dtype=str)
par = par.astype(float)
m = par[0]
b = par[1]

T_data_N = np.loadtxt('T_daten_hoch_mit_uns', dtype=str)
T_data_H = np.loadtxt('T_daten_niedrig_mit_uns', dtype=str)

T_data_N = T_data_N.astype(float)
T_data_H = T_data_H.astype(float)

T_data_N_uns = T_data_N[:,1]
T_data_N = T_data_N[:,0]

T_data_H_uns = T_data_H[:,1]
T_data_H = T_data_H[:,0]

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
    #data_s[i,:] = np.abs(data_N[3*i,:] - data_N[3*i+2,:])/2

data_Cu_N = data_m[:,2] 
#data_Cu_err_N = data_s[:,2]

data_m_ = np.zeros((len(data_H)//3, data_H.shape[1]))
data_s_ = np.zeros((len(data_H)//3, data_H.shape[1]))


for i in range(len(data_H)//3):
    data_m_[i,:] = data_H[3*i+1,:]
    #data_s_[i,:] = np.abs(data_H[3*i,:] - data_H[3*i+2,:])/2
    
data_Cu_H = data_m_[:,2] 
#data_Cu_err_H = data_s_[:,2]

# Bestimmung von beta im nicht-linearen Bereich 

plt.figure(figsize=(20,10))
plt.errorbar(T_data_N_log, np.log(data_Cu_N), xerr= T_data_N_uns / T_data_N ,yerr=np.abs(1/data_Cu_N * sigma_Cu), marker='.', linestyle='', label='data points (pt temp. scale)', zorder=3)
plt.errorbar(T_data_H_log, np.log(data_Cu_H), xerr= T_data_H_uns / T_data_H , yerr=np.abs(1/data_Cu_H * sigma_Cu), marker='.', linestyle='', label='data points (c temp. scale)', zorder=3)
plt.fill_between((T_data_N_log[0]+0.1, T_data_N_log[-2]+0.3) ,-2, 3.2, alpha=0.2, label='linear regime')
plt.fill_between((T_data_H_log[20]+0.4, T_data_H_log[-4]+0.1) ,-2, 3.2, alpha=0.2, label='non-linear regime')
plt.xlabel(r'$\ln(T)$ / K')
plt.ylabel(r'$\ln(R)$ / $\Omega$')
plt.title('Logarithmic plot of the copper resistance as a function of the temperature')
plt.legend()
plt.grid()
plt.show()

#-> nicht linearer Bereich zwischen 3.5 und 4.5 

non_linear_data_N = data_Cu_N[-3:-1]
non_linear_data_H = data_Cu_H[19:22]

non_linear_data = np.append(non_linear_data_N, non_linear_data_H)
non_linear_data = np.sort(non_linear_data)
 
T_data_non_linear_N = T_data_N[-3:-1]
T_data_non_linear_H = T_data_H[19:22]

T_data_non_linear = np.concatenate([T_data_non_linear_N, T_data_non_linear_H])
sort_indices = np.argsort(T_data_non_linear)
T_data_non_linear = T_data_non_linear[sort_indices]

T_data_non_linear_uns_N =  T_data_N_uns[-3:-1]
T_data_non_linear_uns_H =  T_data_H_uns[19:22]

T_data_non_linear_uns = np.concatenate([T_data_non_linear_uns_N, T_data_non_linear_uns_H])

T_data_non_linear_uns = T_data_non_linear_uns[sort_indices]


def non_linear(beta, T_log):
    return beta[0] + beta[1] * T_log

my_model = odr.Model(non_linear)

my_data = odr.RealData(
    x = np.log(T_data_non_linear), 
    y = np.log(non_linear_data),
    sx =  T_data_non_linear_uns / T_data_non_linear , 
    sy = sigma_Cu / non_linear_data    
)       

my_odr = odr.ODR(my_data, my_model, beta0=[1.0, 3])

output = my_odr.run()

sy = sigma_Cu / non_linear_data
sx = T_data_non_linear_uns / T_data_non_linear
steigung = output.beta[1]

plt.figure(figsize=(20,10), layout='tight')
gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 0.25], hspace=0)

plt.subplot(gs[0])
plt.ylabel(r'$\ln(R)$')
plt.title('Logarithmic plot of the non-linear regime (copper)')
plt.grid()
plt.errorbar(np.log(T_data_non_linear), np.log(non_linear_data), yerr = sy, xerr = sx, marker='.', linestyle='', label='data points non-linear regime', zorder=3)
plt.plot(np.log(T_data_non_linear), non_linear(output.beta, np.log(T_data_non_linear)), label='non-linear fit', color='red')
plt.legend()


yerr_res = np.sqrt(sy**2 + (steigung * sx)**2)

plt.subplot(gs[1], sharex=plt.subplot(gs[0]))
residue_non_lin = np.log(non_linear_data) - non_linear(output.beta, np.log(T_data_non_linear))
plt.errorbar(np.log(T_data_non_linear), residue_non_lin, yerr=yerr_res, xerr = sx, fmt='.', label='residuals non-linear fit', zorder=3)
plt.axhline(0, color='black', linestyle='dashed')
plt.xlabel(r'$\ln(T)$')
plt.ylabel('residuals')
plt.grid()
plt.legend()
plt.show()

print("beta = ", output.beta[1])

# -------------------------------------------------------------
# oder Fit an die Kurve: 

def non_linear_2(beta, T):
    return beta[0] * T**beta[1] + beta[2]


my_model = odr.Model(non_linear_2)

my_data = odr.RealData(
    x = T_data_non_linear, 
    y = non_linear_data,
    sx =  T_data_non_linear_uns , 
    sy = sigma_Cu    
)       

my_odr = odr.ODR(my_data, my_model, beta0=[0.1, 6, 0.34])

output_2 = my_odr.run()

sy_2 = sigma_Cu   
sx_2 = T_data_non_linear_uns
beta_1 = output_2.beta[1]
beta_0 = output_2.beta[0]

T_n_l = np.linspace(40, 90, 100)

plt.figure(figsize=(20,10), layout='tight')
gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 0.25], hspace=0)

plt.subplot(gs[0])
plt.ylabel(r'$R$ / $\Omega$')
plt.title('Plot of the non-linear regime (copper)')
plt.grid()
plt.errorbar(T_data_non_linear, non_linear_data, yerr = sy_2, xerr = sx_2, marker='.', linestyle='', label='data points non-linear regime', zorder=3)
plt.plot(T_n_l, non_linear_2(output_2.beta, T_n_l), label='non-linear fit', color='red')
plt.legend()

yerr_res_2 = np.sqrt(sy_2**2 + (beta_1 * beta_0 * T_data_non_linear **(beta_1-1) * sx_2)**2)

plt.subplot(gs[1], sharex=plt.subplot(gs[0]))
residue_non_lin_2 = non_linear_data - non_linear_2(output_2.beta, T_data_non_linear)
plt.errorbar(T_data_non_linear, residue_non_lin_2, yerr=yerr_res_2, xerr = sx_2, fmt='.', label='residuals non-linear fit', zorder=3)
plt.axhline(0, color='black', linestyle='dashed')
plt.xlabel(r'$T$ / K')
plt.ylabel('residuals')
plt.grid()
plt.legend()
plt.show()

print("beta = ", output_2.beta[1])


# Bestimmmung des linearen Widerstandskoeffizienten alpha 
# durch lineare Regression der Daten im linearen Bereich

linear_data_N = data_Cu_N[:-3]
linear_data_H = data_Cu_H[:19]

linear_data_T_N = T_data_N[:-3]
linear_data_T_H = T_data_H[:19]

linear_data_T_N_uns = T_data_N_uns[:-3]
linear_data_T_H_uns = T_data_H_uns[:19]


#beta[0] = R_0, beta[1] = alpha

def linear(beta, T):
    return beta[0]* (beta[1] * T + 1)

my_model = odr.Model(linear)

my_data = odr.RealData(
    x = linear_data_T_N-273.15, 
    y = linear_data_N,
    sx = linear_data_T_N_uns, 
    sy = sigma_Cu    
)       

my_odr = odr.ODR(my_data, my_model, beta0=[20, 17/20])

output_linear = my_odr.run()

sx = linear_data_T_N_uns
sy = sigma_Cu
beta_1 = output_linear.beta[1]
beta_0 = output_linear.beta[0]

my_data_H = odr.RealData(
    x = linear_data_T_H-273.15, 
    y = linear_data_H,
    sx = linear_data_T_H_uns, 
    sy = sigma_Cu    
)

my_odr_H = odr.ODR(my_data_H, my_model, beta0=[20, 17/20])

output_linear_H = my_odr_H.run()
sx_H = linear_data_T_H_uns
sy_H = sigma_Cu
beta_1_H = output_linear_H.beta[1]
beta_0_H = output_linear_H.beta[0]
print(output_linear.beta, output_linear_H.beta)
plt.figure(figsize=(20,10), layout='tight')

gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 0.25], hspace=0)

plt.subplot(gs[0])
plt.errorbar(linear_data_T_N-273.15, linear_data_N, yerr=sy, xerr= sx, marker='.', linestyle='', label='data points (pt temp. scale)', zorder=3)
plt.errorbar(linear_data_T_H-273.15, linear_data_H, yerr=sy, xerr= sx_H, marker='.', linestyle='', label='data points (c temp. scale)', zorder=3)

T_ = np.linspace(-175, 20, 1000)
plt.plot(T_, linear(output_linear.beta, T_), label='linear fit (pt temp. scale)', color='blue')
plt.plot(T_, linear(output_linear_H.beta, T_), label='linear fit (c temp. scale)', color='orange')
plt.ylabel(r'$R$ / $\Omega$')
plt.title('Plot of the copper resistance as a function of the temperature')
plt.grid()
plt.legend()    

res_yerr_N = np.sqrt(sy**2 + (beta_1 * beta_0  * sx)**2)
res_yerr_H = np.sqrt(sy_H**2 + (beta_1_H * beta_0_H * sx_H)**2)


plt.subplot(gs[1], sharex=plt.subplot(gs[0]))
residue_N = linear_data_N - linear(output_linear.beta, linear_data_T_N-273.15)
residue_H = linear_data_H - linear(output_linear_H.beta, linear_data_T_H-273.15)

plt.errorbar(linear_data_T_N-273.15, residue_N, yerr=res_yerr_N, xerr= sx, fmt='.', label='residuals (pt temp. scale)', zorder=3)
plt.errorbar(linear_data_T_H-273.15, residue_H, yerr=res_yerr_H, xerr=sx_H, fmt='.', label='residuals (c temp. scale)', zorder=3)
plt.axhline(0, color='black', linestyle='dashed')
plt.xlabel(r'$T$ / °C')
plt.ylabel('residuals')
plt.grid()
plt.legend()
plt.show()


print("linearer Widerstandskoeffizient Fit Pt alpha = ", output_linear.beta[1])
print("Fit C alpha = ", output_linear_H.beta[1], "\n alpha aus Literatur: ", 0.0039)
#print("R_0 = ", output_linear.beta[0], "\n R_0 aus Fit: ", output_linear_H.beta[0])

#Quelle alpha: https://www.chemie.de/lexikon/Temperaturkoeffizient.html#:~:text=Table_title:%20Temperaturkoeffizient%20des%20elektrischen%20Widerstands%20Table_content:%20header:,%C2%B7%2010%2D3%20%7C%20%CE%B1%20in%20K%2D1:%20%7C




'''
# Fit mit Gerade
popt_non_lin, cov_non_lin = curve_fit(non_linear, np.log(T_data_non_linear), np.log(non_linear_data), p0=[-1, 0.1])

plt.figure(figsize=(20,10), layout='tight')

gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 0.25], hspace=0)

plt.subplot(gs[0])
plt.ylabel(r'$\ln(R)$')
plt.title('Logarithmic plot of the non-linear regime (copper)')
plt.legend()
plt.grid()
plt.errorbar(np.log(T_data_non_linear), np.log(non_linear_data), yerr=sy, xerr=sx, marker='.', linestyle='', label='data points non-linear regime', zorder=3)
plt.plot(np.log(T_data_non_linear), non_linear(np.log(T_data_non_linear), *popt_non_lin), label='non-linear fit')

plt.subplot(gs[1], sharex=plt.subplot(gs[0]))
residue_non_lin = np.log(non_linear_data) - non_linear(np.log(T_data_non_linear), *popt_non_lin)
plt.errorbar(np.log(T_data_non_linear), residue_non_lin, yerr=yerr_res, xerr=sx, fmt='.', label='residuals non-linear fit', zorder=3)
plt.axhline(0, color='black', linestyle='dashed')
plt.xlabel(r'$\ln(T)$')
plt.ylabel('residuals')
plt.grid()
plt.legend()
plt.show()

#Fit mit Kurve 

T_n_l = np.linspace(40, 90, 100)

popt_non_lin_2, cov_non_lin_2 = curve_fit(non_linear_2, T_data_non_linear, non_linear_data, p0=[1, 3, 1])

plt.figure(figsize=(20,10), layout='tight')

gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 0.25], hspace=0)

plt.subplot(gs[0])
plt.errorbar(T_data_non_linear, non_linear_data, yerr=np.abs(sigma_Cu), marker='.', linestyle='', label='data points non-linear regime', zorder=3)
plt.plot(T_n_l, non_linear_2(T_n_l, *popt_non_lin_2), label='non-linear fit')
plt.ylabel(r'$R$ / $\Omega$')
plt.title('Plot of the non-linear regime (copper)')
plt.legend()
plt.grid()


plt.subplot(gs[1], sharex=plt.subplot(gs[0]))
residue_non_lin_2 = non_linear_data - non_linear_2(T_data_non_linear, *popt_non_lin_2)
plt.errorbar(T_data_non_linear, residue_non_lin_2, yerr=np.abs(sigma_Cu), fmt='.', label='residuals non-linear fit', zorder=3)
plt.axhline(0, color='black', linestyle='dashed')
plt.xlabel(r'$T$ / K')
plt.ylabel('residuals')
plt.grid()
plt.legend()
plt.show()
print("beta = ", popt_non_lin_2[1])



popt_N, cov_N = curve_fit(linear, linear_data_T_N-273.15, linear_data_N, p0=[100, 0.0039])
#popt_H, cov_H = curve_fit(linear, T_data_H[:-5], data_Cu_H[:-5])

plt.figure(figsize=(20,10), layout='tight')

gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 0.25], hspace=0)

plt.subplot(gs[0])
plt.errorbar(linear_data_T_N-273.15, linear_data_N, yerr=np.abs(sigma_Cu), marker='.', linestyle='', label='data points (pt temp. scale)', zorder=3)
#plt.errorbar(T_data_H-273.15, data_Cu_H, yerr=np.abs(sigma_Cu), marker='.', linestyle='', label='data points (c temp. scale)', zorder=3)

T_ = np.linspace(-160, 20, 1000)

plt.plot(T_, linear(T_, *popt_N), label='linear fit (pt temp. scale)')
#plt.plot(T_data_H-273.15, linear_H(T_data_H-273.15, *popt_H), label='linear fit (c temp. scale)')
plt.ylabel(r'$R$ / $\Omega$')
plt.title('Plot of the copper resistance as a function of the temperature')
plt.legend()
plt.grid()

'''