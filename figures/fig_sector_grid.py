"""
fig_sector_grid.py
------------------
Figure: The four-sector carry grid for K = 4.

For every pair (X, Y) of 4-bit numbers (8..15), determines:
  - Sector: (high_X, high_Y) where high = 1 iff second-highest bit = 1
  - D-odd status: product X*Y < 2^(2K-1) = 128
  - Carry weight via the cascade valuation

The grid is 8x8 (X = 8..15, Y = 8..15), coloured by:
  - D-even pairs: light grey
  - Sector (1,1): dark grey (always D-even — Theorem 4)
  - D-odd pairs in (0,0), (0,1), (1,0): warm/cool colour by carry weight (-1..+max)
"""

import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

K = 4
N = 1 << K          # 16
HALF = N >> 1       # 8  — numbers in [8, 15]
D_ODD_MAX = 1 << (2 * K - 1)   # 128  — products below this are D-odd


def carry_chain(x, y, K):
    """Return the carry chain c[0..2K] for x*y (K-bit x K-bit)."""
    bits_x = [(x >> i) & 1 for i in range(K)]
    bits_y = [(y >> i) & 1 for i in range(K)]
    # column sums s[j] = sum_{i} bits_x[i] * bits_y[j-i]  for valid i
    s = [0] * (2 * K)
    for i in range(K):
        for jj in range(K):
            s[i + jj] += bits_x[i] * bits_y[jj]
    # propagate carries
    c = [0] * (2 * K + 1)
    for j in range(2 * K):
        total = s[j] + c[j]
        c[j + 1] = total >> 1
    return c


def cascade_weight(c, K):
    """
    Cascade valuation: find the highest position M where c[M] > 0,
    then read c[M-1]; weight = c[M-1] - 1.
    Returns None for D-even pairs (no bridge).
    """
    # D-odd: c starts and ends at 0 with a hump
    # highest nonzero carry position
    M = None
    for j in range(2 * K, -1, -1):
        if c[j] > 0:
            M = j
            break
    if M is None or M == 0:
        return None
    val = c[M - 1] - 1   # weight: 0 → -1, 1 → 0, 2 → +1, …
    return val


# ── compute grid ──────────────────────────────────────────────────────────────
grid_weight = np.full((8, 8), np.nan)   # NaN = D-even
grid_sector = np.zeros((8, 8), dtype=int)  # 0=00,1=01,2=10,3=11

for ix, x in enumerate(range(HALF, N)):
    for iy, y in enumerate(range(HALF, N)):
        hx = (x >> (K - 2)) & 1   # second-highest bit
        hy = (y >> (K - 2)) & 1
        sector = hx * 2 + hy       # 0,1,2,3
        grid_sector[ix, iy] = sector
        prod = x * y
        if prod < D_ODD_MAX:       # D-odd
            c = carry_chain(x, y, K)
            w = cascade_weight(c, K)
            if w is not None:
                grid_weight[ix, iy] = w

# ── plotting ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 6))

# sector background colours (light)
sector_bg = {
    0: "#dbeafe",   # 00: light blue
    1: "#dcfce7",   # 01: light green
    2: "#fef9c3",   # 10: light yellow
    3: "#e5e7eb",   # 11: grey (always D-even)
}
for ix in range(8):
    for iy in range(8):
        s = grid_sector[ix, iy]
        rect = plt.Rectangle([iy - 0.5, ix - 0.5], 1, 1,
                              facecolor=sector_bg[s], edgecolor="white",
                              linewidth=0.8, zorder=1)
        ax.add_patch(rect)

# D-odd carry weight as filled circles
w_min = -1
w_max = int(np.nanmax(grid_weight)) if not np.all(np.isnan(grid_weight)) else 1
cmap = plt.get_cmap("RdYlBu_r")

for ix in range(8):
    for iy in range(8):
        w = grid_weight[ix, iy]
        if not np.isnan(w):
            norm_w = (w - w_min) / max(w_max - w_min, 1)
            color = cmap(norm_w)
            circle = plt.Circle((iy, ix), 0.32, color=color, zorder=3)
            ax.add_patch(circle)
            ax.text(iy, ix, str(int(w)), ha="center", va="center",
                    fontsize=7, fontweight="bold", color="white" if norm_w < 0.3 or norm_w > 0.7 else "black",
                    zorder=4)
        elif grid_sector[ix, iy] == 3:
            # (1,1) sector — explicitly mark as theorem 4
            ax.text(iy, ix, "✕", ha="center", va="center",
                    fontsize=10, color="#9ca3af", zorder=4)
        else:
            # D-even in non-(1,1) sector
            ax.text(iy, ix, "·", ha="center", va="center",
                    fontsize=14, color="#9ca3af", zorder=4)

# sector labels (outside grid)
sector_labels = {0: "Sector (0,0)", 1: "Sector (0,1)", 2: "Sector (1,0)", 3: "Sector (1,1)"}
# draw sector boundary lines
ax.axvline(1.5, color="#374151", lw=1.5, zorder=5)
ax.axhline(1.5, color="#374151", lw=1.5, zorder=5)

# axis labels
x_vals = list(range(HALF, N))
y_vals = list(range(HALF, N))
ax.set_xticks(range(8))
ax.set_xticklabels([str(v) for v in y_vals], fontsize=9)
ax.set_yticks(range(8))
ax.set_yticklabels([str(v) for v in x_vals], fontsize=9)
ax.set_xlabel("$Y$ (second factor)", fontsize=11)
ax.set_ylabel("$X$ (first factor)", fontsize=11)
ax.set_title(f"Carry weight grid for $K = {K}$ (4-bit × 4-bit)\nCircles = D-odd pairs, colour = cascade carry weight", fontsize=11)
ax.set_xlim(-0.5, 7.5)
ax.set_ylim(-0.5, 7.5)
ax.set_aspect("equal")

# legend
patches = [
    mpatches.Patch(facecolor="#dbeafe", edgecolor="#374151", label="Sector (0,0)  $X$ low, $Y$ low"),
    mpatches.Patch(facecolor="#dcfce7", edgecolor="#374151", label="Sector (0,1)  $X$ low, $Y$ high"),
    mpatches.Patch(facecolor="#fef9c3", edgecolor="#374151", label="Sector (1,0)  $X$ high, $Y$ low"),
    mpatches.Patch(facecolor="#e5e7eb", edgecolor="#374151", label="Sector (1,1)  always D-even  (Thm 4)"),
]
from matplotlib.lines import Line2D
legend_circ = Line2D([0], [0], marker="o", color="w", markerfacecolor=cmap(0.0),
                     markersize=10, label="weight $= -1$")
legend_circ2 = Line2D([0], [0], marker="o", color="w", markerfacecolor=cmap(0.5),
                      markersize=10, label="weight $= 0$")
legend_circ3 = Line2D([0], [0], marker="o", color="w", markerfacecolor=cmap(1.0),
                      markersize=10, label="weight $= +1$")
ax.legend(handles=patches + [legend_circ, legend_circ2, legend_circ3],
          fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)

plt.tight_layout()
out = "papers/carry-arithmetic-P1-pi-spectral/figures/fig_sector_grid.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")
