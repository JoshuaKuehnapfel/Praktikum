#%%

import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy as unp

P0 = ufloat(23.19, 0.30)
PP = ufloat(0.831, 0.008)
PB = ufloat(0.077, 0.005)

Verlust_PP = 1-PP/P0
Verlust_PB = 1-PB/P0
print(f"Verlust PP: {Verlust_PP*100:.3f} %")
print(f"Verlust PB: {Verlust_PB*100:.3f} %")


