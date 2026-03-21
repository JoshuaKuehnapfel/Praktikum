import numpy as np
import os
os.chdir("Praktikum/GEM/Daten")

name = ['Kadmium', 'Eisen', 'xray']
data_rausch = np.loadtxt('Rauschspectrum.csv', delimiter=',', skiprows=1)

for n in name:
    data = np.loadtxt(n + 'spectrum.csv', delimiter=',', skiprows=1)-data_rausch
    print(n)
    print('Summe: ', np.sum(data).round(2))