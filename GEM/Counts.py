<<<<<<< HEAD
import numpy as np
import os
os.chdir("Praktikum/GEM/Daten")

name = ['Kadmium', 'Eisen', 'xray']
data_rausch = np.loadtxt('Rauschspectrum.csv', delimiter=',', skiprows=1)

for n in name:
    data = np.loadtxt(n + 'spectrum.csv', delimiter=',', skiprows=1)-data_rausch
    print(n)
    print('Summe: ', np.sum(data).round(2))
=======
import numpy as np

daten = np.loadtxt('GEM/Daten/Eisenspektrum.csv', skiprows = 1, delimiter = ',')
summe = np.sum(daten[1,:])
print('Die Summe der Zählraten beträgt: ', summe)
>>>>>>> 4459621dedf4e1a401f57a4774b5e6cb8f78cd3b
