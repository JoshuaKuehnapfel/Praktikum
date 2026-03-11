import numpy as np

daten = np.loadtxt('GEM/Daten/Eisenspektrum.csv', skiprows = 1, delimiter = ',')
summe = np.sum(daten[1,:])
print('Die Summe der Zählraten beträgt: ', summe)
