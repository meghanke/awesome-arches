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


# 1)  Load photometry

CSV_PATH = Path(__file__).with_name("arches_cluster.csv")
df = pd.read_csv(CSV_PATH)


# 2)  Compute color index and filter valid rows

color_index = df["mF127M (mag)"] - df["mF153M (mag)"]
mask        = color_index.notna() & df["mF153M (mag)"].notna()
df          = df[mask].copy()
color_index = color_index[mask]
mag153      = df["mF153M (mag)"]


# 3)  Marker sizes: larger points for brighter (smaller-mag) stars
#     Simple scaling:  size ∝ 10^(-0.4 Δm)

base_size = 1500                              
m_min     = mag153.min()
sizes     = base_size * 10 ** (-0.4 * (mag153 - m_min))

# 4)  Plot

fig, ax = plt.subplots(figsize=(7, 8))

sc = ax.scatter(
    color_index,          # x
    mag153,               # y
    c=color_index,        # point colors
    cmap="turbo",         # rainbow-like colormap
    s=sizes,              # marker areas
    edgecolors="k",
    linewidths=0.3
)

ax.invert_yaxis()         # bright stars appear at the top

ax.set_xlabel(r"$(m_{127} - m_{153})$  color  (mag)")
ax.set_ylabel(r"$m_{153}$  apparent magnitude (mag)")
ax.set_title("Arches Cluster – Near-IR Color–Magnitude Diagram")
ax.grid(alpha=0.3)

cbar = fig.colorbar(sc, ax=ax, pad=0.02)
cbar.set_label(r"$(m_{127} - m_{153})$  color  (mag)")

plt.tight_layout()
plt.show()