import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit

os.chdir('C:/VS_Studio_Latex/Uni (lLokale Kopie)/Praktikum/T01/Alpha in Luft')
daten_spek = np.loadtxt('Alpha_Luft_Spektrum.TKA', skiprows=2)
daten_peakl = np.loadtxt('Alpha_Luft_2.Peak_links.TKA', skiprows=2)
daten_peakm = np.loadtxt('Alpha_Luft_2.Peak_mitte.TKA', skiprows=2)
daten_peakr = np.loadtxt('Alpha_Luft_2.Peak_rechts.TKA', skiprows=2)

daten = [daten_spek, daten_peakl, daten_peakm, daten_peakr]

def gaussian(x, A, sigma, x0):
    return A * np.exp(-0.5 * ((x - x0) / sigma)**2)

def peak_fitter(data):
    peak_1 = np.argmax(data)    
    data_cleaned = np.delete(data, np.arange(max(0, peak_1 - 500), min(len(data), peak_1 + 500)))
    peak_2 = np.argmax(data_cleaned)
    data_cleaned2 = np.delete(data_cleaned, np.arange(max(0, peak_2 - 1000), min(len(data_cleaned), peak_2 + 1000)))
    return [peak_1, np.where(data == data_cleaned[peak_2])[0][0], np.where(data == np.max(data_cleaned2))[0][0]]    






for i, d in enumerate(daten):
    peaks = peak_fitter(d)
    print(peaks)
    x = np.arange(len(d))
    plt.figure(figsize=(10, 6))
    plt.plot(x, d, label='Daten')    
    plt.vlines(peaks, min(d), max(d), colors='green', linestyles='dashed', label='Gefundene Peaks')
    plt.legend()
    plt.grid()
    plt.show()



lim1 = 7700
lims = 8000
