#%%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks_cwt

daten_spek = np.loadtxt('Alpha_Luft_Spektrum.TKA', skiprows=2)
daten_peakl = np.loadtxt('Alpha_Luft_2.Peak_links.TKA', skiprows=2)
daten_peakm = np.loadtxt('Alpha_Luft_2.Peak_mitte.TKA', skiprows=2)
daten_peakr = np.loadtxt('Alpha_Luft_2.Peak_rechts.TKA', skiprows=2)

daten = [daten_spek, daten_peakl, daten_peakm, daten_peakr]

def gaussian(x, A, sigma, x0):
    return A * np.exp(-0.5 * ((x - x0) / sigma)**2)

peaks = []
    
    




for i, d in enumerate(daten):
    plt.figure(figsize=(10, 6))
    peaks = peak_fitter(d)
    x = np.arange(len(d))
    plt.plot(x, d, label='Daten', alpha = 1)    
    plt.vlines(peaks, min(d), max(d), colors='green', alpha =1, linestyles='dashed', label='Gefundene Peaks')
    plt.legend()
    plt.grid()
    plt.show()



lim1 = 7700
lims = 8000

# %%
