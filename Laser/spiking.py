#%%

from Schnell_auswertung import linear, quadratic, fit_mit_plot
import numpy as np
import matplotlib.pyplot as plt
plt.clf()
Data1 = np.loadtxt('ALL0001/F0001CH1.csv', delimiter=',', skiprows=18, usecols=(3, 4))
Data2 = np.loadtxt('ALL0002/F0002CH1.csv', delimiter=',', skiprows=18, usecols=(3, 4))
Data1[:, 1] += 6.400000e-03
plt.plot(Data1[:, 0], Data1[:, 1], label='Data1')
plt.plot(Data2[:, 0], Data2[:, 1], label='Data2', linewidth=0.5)

Gezoomte_Daten1 = np.where((Data1[:, 0] > 0-6.4e-3) & (Data1[:, 0] < 0.00335))
plt.show()
plt.plot(Data1[Gezoomte_Daten1][:, 0]*1e3, Data1[Gezoomte_Daten1][:, 1], label='Spike-Daten', linewidth=0.5)
plt.hlines(Data1[Gezoomte_Daten1][:, 1][-10], 0, 0.00335*1e3, colors='red', linestyles='dashed', label='CW-Plateau', alpha = 0.5)
plt.hlines(Data1[Gezoomte_Daten1][:, 1].max(), 0, 0.00335*1e3, colors='green', linestyles='dashed', label='Spiking-Peak', alpha = 0.5)
plt.xlabel(r"Zeit [$\mu$s]")
plt.ylabel("Amplitude [V]")
delta = Data1[Gezoomte_Daten1][:, 1].max() - Data1[Gezoomte_Daten1][:, 1][-10]
plt.text(2.5, Data1[Gezoomte_Daten1][:, 1][-10]+delta/2, f"Delta = {delta:.3f} V", fontsize=10, ha='center', va='center', bbox=dict(facecolor='white', alpha=0.5))
plt.annotate('', xy=(2, Data1[Gezoomte_Daten1][:, 1].max()), xytext=(2, Data1[Gezoomte_Daten1][:, 1][-10]), arrowprops=dict(arrowstyle='<->', color='black', lw=0.8))
plt.grid()
plt.legend()
plt.savefig('Spiking.svg', dpi=300)
plt.title("Spiking-Effekt beim Nd:YAG-Laser")
plt.show()
delta2 = Data1[Gezoomte_Daten1][:, 1][-10] - Data1[Gezoomte_Daten1][:, 1].min()
print(f"Delta2 = {delta2:.3f} V")