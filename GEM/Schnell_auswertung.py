import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scopt
import uncertainties as unc
import scipy.constants as scon

U = [3100, 3400, 3300, 3200, 3350, 3000, 3250]
U0 = [2900, 3000, 3100, 3200, 3300, 3350, 3400]
U_80_40 = np.array([2800, 2880, 2960, 3040, 3120, 3210, 3300])
R_an = [2207, 14639, 13502, 7456, 14360, 2, 12109]
R_aus = [3, 39, 20, 5, 30, 0, 18]
x,y = np.loadtxt('V_Drift_80_20.csv', delimiter=',',  unpack=True, usecols=(1,2))

x0,y0 = np.loadtxt('V_Drift.csv', delimiter=',',  unpack=True, usecols=(1,2))
fig, ax = plt.subplots(2, 2)

N = unc.ufloat(322, 2.8)
e = scon.elementary_charge
R = [unc.ufloat(i, j) for i, j in zip(R_an, R_aus)]


ax[0,0].plot(U_80_40, -x[::2], 'x', label='Messwerte')
ax[0,0].set_xlabel('Spannung U [V]')
ax[0,0].set_ylabel('I_RO [A]')
ax[0,0].set_title('80/20')
ax[0,0].legend()


Gain = [R_i / N / e for R_i in R]
print([g.n for g in Gain])
ax[0,1].plot(U, [g.n for g in Gain], 'x', label='Verstärkung')
ax[0,1].set_xlabel('Spannung U [V]')
"""ax[1,1].plot(U, R_an.n, 'x', label='Zählung (Quelle an)')
ax[1,1].plot(U, R_aus.n, 'x', label='Zählung (Quelle aus)')"""
ax[1,1].set_xlabel('Spannung U [V]')
ax[1,1].set_ylabel('Rate [Hz]')
ax[1,1].legend()

ax[1,0].plot(U0, x0[::2], 'x', label='Messwerte')
ax[1,0].set_xlabel('Spannung U [V]')
ax[1,0].set_ylabel('I_RO [A]')
ax[1,0].set_title('70/30')
plt.tight_layout()
plt.show()
