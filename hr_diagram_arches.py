# hr_diagram_arches.py
#
# H–R–style color–magnitude diagram for the Arches cluster
#  • x-axis  : (mF127M – mF153M) ≈ J–H color
#  • y-axis  :  mF153M apparent magnitude  (brighter = up)
#  • marker  : point size ∝ brightness; color follows a rainbow gradient

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path        
# ------------------------------------------------------------
# 1)   user-supplied constants
# ------------------------------------------------------------
DIST_PC = 8000            # distance to the Arches cluster (≈ 8 kpc)
BC_153  = -1.5            # crude bolometric correction for F153M
Mbol_sun = 4.74
T_sun    = 5778.0         # K

# ------------------------------------------------------------
# helper: distance modulus → absolute magnitude
# ------------------------------------------------------------
def apparent2absolute(m, d_pc):
    """m → M given distance in parsec."""
    return m - 5*np.log10(d_pc/10)

# ------------------------------------------------------------
# helper: Stefan-Boltzmann –> radius (in solar units)
# ------------------------------------------------------------
def radius_from_L_T(L_over_Lsun, Teff):
    """Return R/Rsun from L/Lsun and T."""
    return np.sqrt(L_over_Lsun) / (Teff / T_sun)**2

# -----------------------------------------------------------------
# 2)  load and pre-process exactly as you did
# -----------------------------------------------------------------
CSV_PATH = Path(__file__).with_name("arches_cluster.csv")
df = pd.read_csv(CSV_PATH)

color_index = df["mF127M (mag)"] - df["mF153M (mag)"]
mask        = color_index.notna() & df["mF153M (mag)"].notna()
df          = df[mask].copy()
color_index = color_index[mask]
m153        = df["mF153M (mag)"]

# -----------------------------------------------------------------
# 3)  Apparent  CMD  
# -----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 7))

sizes = 1500 * 10**(-0.4 * (m153 - m153.min()))
sc = ax.scatter(color_index, m153, c=color_index, cmap="turbo",
                s=sizes, edgecolors="k", linewidths=0.3)

ax.invert_yaxis()
ax.set_xlabel(r"$(m_{127}-m_{153})$ (mag)")
ax.set_ylabel(r"$m_{153}$ (mag)")
ax.set_title("Arches cluster – apparent CMD")
fig.colorbar(sc, ax=ax, pad=0.02, label="colour (mag)")

# -----------------------------------------------------------------
# 4)  Absolute CMD  (apply distance modulus)
# -----------------------------------------------------------------
M153 = apparent2absolute(m153, DIST_PC)

fig2, ax2 = plt.subplots(figsize=(6, 7))
sc2 = ax2.scatter(color_index, M153, c=color_index, cmap="turbo",
                  s=sizes, edgecolors="k", linewidths=0.3)

ax2.invert_yaxis()
ax2.set_xlabel(r"$(m_{127}-m_{153})$ (mag)")
ax2.set_ylabel(r"$M_{153}$ (mag)")
ax2.set_title(f"Absolute CMD  (d = {DIST_PC/1000:.1f} kpc)")
fig2.colorbar(sc2, ax=ax2, pad=0.02, label="colour (mag)")

# -----------------------------------------------------------------
# 5)  *Very* rough temperature & luminosity estimates -------------
#      (good enough for demonstration)
# -----------------------------------------------------------------
# crude colour→Teff calibration: Teff ≈ a – b*(colour)
a, b = 13000, 8500         # picked to map colour~0 → 13 kK, colour~2 → ~1 kK
Teff = a - b*color_index.clip(0, 1.5)   # keep within 0–1.5 mag
# bolometric magnitude
Mbol = M153 + BC_153
# L/Lsun from Mbol
L_over_Lsun = 10**((Mbol_sun - Mbol)/2.5)
# radius via Stefan-Boltzmann
R_over_Rsun = radius_from_L_T(L_over_Lsun, Teff)

# -----------------------------------------------------------------
# 6)  H-R diagram with constant-R lines
# -----------------------------------------------------------------
fig3, ax3 = plt.subplots(figsize=(7, 6))
logT = np.log10(Teff)
logL = np.log10(L_over_Lsun)

sc3 = ax3.scatter(logT, logL, c=R_over_Rsun, cmap="inferno",
                  s=30, edgecolors="none")
ax3.set_xlim(ax3.get_xlim()[::-1])  # hottest (logT big) on left
ax3.set_xlabel(r"$\log\,T_{\rm eff}$  [K]")
ax3.set_ylabel(r"$\log\,L/L_\odot$")
ax3.set_title("H-R diagram – Arches cluster")

# constant-radius loci
for R in [0.5, 1, 5, 20, 100]:
    logT_grid = np.linspace(3.6, 4.6, 200)
    logL_grid = 2*np.log10(R) + 4*(logT_grid - np.log10(T_sun))
    ax3.plot(logT_grid, logL_grid, ls="--",
             label=f"{R:g} $R_\\odot$", lw=1)

ax3.legend(title="constant R")
cbar3 = fig3.colorbar(sc3, ax=ax3, pad=0.03)
cbar3.set_label(r"estimated $R/R_\odot$")

plt.tight_layout()
plt.show()