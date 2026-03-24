#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties import correlated_values 
from scipy.special import erf
x1 = 9.05-5
x2 = 9.05-4.5
def linear(x, a, b):
    return a*x + b
y1 = 0.115
y1_unc = 0.010
y2 = 0.0125
y2_unc = 0.0015

popt, pcov = curve_fit(linear, [x1, x2], [y1, y2], sigma=[y1_unc, y2_unc], absolute_sigma=True)
a, b = popt 
a_unc, b_unc = np.sqrt(np.diag(pcov))
a = ufloat(a, a_unc)
b = ufloat(b, b_unc)
# a und b als korrelierte Unsicherheitsgrößen aus Fit-Ergebnis erzeugen
a_corr, b_corr = correlated_values(popt, pcov)

x_0 = -b_corr / a_corr
print(x_0)
R_lit = ufloat(7.06, 0.35)
print(R_lit - x_0)

x3 = 9.05-6
x4 = 9.05-6.5
y3 = 0.360
y3_unc = 0.015
y4 = 0.365
y4_unc = 0.015

x_data = np.array([x1, x2, x3, x4], dtype=float)
y_data = np.array([y1, y2, y3, y4], dtype=float)
y_sigma = np.array([y1_unc, y2_unc, y3_unc, y4_unc], dtype=float)

def err_func(x, A, x0, s, C):
    return C + 0.5 * A * (1.0 + erf((x - x0) / (np.sqrt(2.0) * s)))

p0 = [1, 3.5, 1, -0.2]
# Startwerte passend zu den aktuell fallenden Daten
C0 = float(np.max(y_data))                 # oberes Plateau
A0 = float(np.min(y_data) - C0)            # negativer Sprung
y_half = C0 + 0.5 * A0
x00 = float(x_data[np.argmin(np.abs(y_data - y_half))])
s0 = max((x_data.max() - x_data.min()) / 8.0, 0.05)

p0 = [A0, x00, s0, C0]
bounds = ([-np.inf, x_data.min() - 2.0, 1e-6, -np.inf],
          [ np.inf, x_data.max() + 2.0, np.inf,  np.inf])
bounds = ([-np.inf, x_data.min() - 2.0, 1e-6, -np.inf],
          [0, x_data.max() + 2.0, np.inf, np.inf])

popt_erf, pcov_erf = curve_fit(
    err_func, x_data, y_data, p0=p0, sigma=y_sigma, absolute_sigma=True, bounds=bounds
)

A, x0_erf, s_erf, C = correlated_values(popt_erf, pcov_erf)
print(f"A = {A}")
print(f"x0 = {x0_erf}")
print(f"s  = {s_erf}")
print(f"C  = {C}")

x_fit = np.linspace(x_data.min() - 0.3, x_data.max() + 0.3, 400)
plt.errorbar(x_data, y_data, yerr=y_sigma, fmt="o", capsize=3, label="data")
plt.plot(x_fit, err_func(x_fit, *popt_erf), label="erf fit")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.tight_layout()
plt.show()
# Residuenplot (zusammen mit Fit)
y_fit_data = err_func(x_data, *popt_erf)
residuals = y_data - y_fit_data

fig, (ax_top, ax_res) = plt.subplots(
    2, 1, sharex=True, figsize=(6, 6),
    gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05}
)

# Oberes Panel: Daten + Fit

ax_top.grid(True)
ax_res.grid(True)
ax_res.errorbar(x_data, residuals, yerr=y_sigma, fmt="o", capsize=3)
ax_res.axhline(0, color="black", lw=1.5, ls="--")
ax_res.set_xlabel("Abstand [cm]")
ax_res.set_ylabel("residuen")
plt.tight_layout()



x_wende = popt_erf[1]
x_wende_u = x0_erf  

ax_top.errorbar(x_data, y_data, yerr=y_sigma, fmt="o", capsize=3, label="Messdaten mit Unsicherheiten")
ax_top.plot(x_fit, err_func(x_fit, *popt_erf), color="C1", label="erf-Fit")

ax_top.axvline(x_wende, color="crimson", ls="--", lw=1.8, label="Wendepunkt")
x_w_min = x_wende_u.n - x_wende_u.s
x_w_max = x_wende_u.n + x_wende_u.s

y_min = np.min(y_data - y_sigma) - 0.02
y_max = np.max(y_data + y_sigma) + 0.02

ax_top.fill_between(
    x_fit, y_min, y_max,
    where=(x_fit >= x_w_min) & (x_fit <= x_w_max),
    color="crimson", alpha=0.18, interpolate=True,
    label=r"Unsicherheit $\overline{R} \pm R_\sigma$"
)
textbox = f"mittlere Reichweite\n $\\overline{{R}} = {x_wende_u.n:.3f} \\pm {x_wende_u.s:.3f}$" + ' cm'
ax_top.text(
    0.7, 0.4, textbox,
    transform=ax_top.transAxes,
    fontsize=10,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="gray")
)


ax_top.set_xlabel("Abstand [cm]")
ax_top.set_ylabel("Strom [0.1 nA]")
ax_top.grid(True, alpha=0.3)
ax_top.legend()
fig.suptitle(r"Bestimmung der mittlere Reichweite der $^{214}$Po-$\alpha$-Teilchen", fontsize=14)
plt.tight_layout()
plt.savefig("Szintilator_fit.svg", dpi=300)