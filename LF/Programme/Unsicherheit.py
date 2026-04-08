#%% 
import numpy as np 

data = np.loadtxt('Rauschmessung', dtype=str)

data = np.char.replace(data, ',', '.')

data = data.astype(float)

print(np.std(data[:,0], ddof=1))

