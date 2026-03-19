#!/usr/bin/env python3
"""
fig_sector_layout.py — Visual layout of the four sectors for K=4 binary multiplication.

Shows the 2×2 sector grid defined by the second-highest bit of X and Y,
with D-odd pair counts and the (1,1) exclusion (Theorem 4).

Output: fig_sector_layout.png
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

fig, ax = plt.subplots(figsize=(7, 6))

# Grid dimensions
x0, y0 = 0.1, 0.1
w, h = 0.4, 0.4

# Sector data: (row, col) -> (label, d_odd_count, color, symbol)
sectors = {
    (0, 0): ("Sector (0,0)", "16 D-odd pairs", "#4A90D9", "✓"),
    (0, 1): ("Sector (0,1)", "8 D-odd pairs",  "#7BB3E0", "✓"),
    (1, 0): ("Sector (1,0)", "8 D-odd pairs",  "#E8A838", "✓"),
    (1, 1): ("Sector (1,1)", "0 D-odd pairs",  "#D0D0D0", "✗"),
}

for (r, c), (label, count, color, sym) in sectors.items():
    xpos = x0 + c * w
    ypos = y0 + (1 - r) * h  # row 0 on top

    rect = patches.FancyBboxPatch(
        (xpos + 0.01, ypos + 0.01), w - 0.02, h - 0.02,
        boxstyle="round,pad=0.02",
        facecolor=color, edgecolor="#333333", linewidth=2,
        alpha=0.85
    )
    ax.add_patch(rect)

    # Sector label
    ax.text(xpos + w/2, ypos + h * 0.65, label,
            ha='center', va='center', fontsize=13, fontweight='bold',
            color='white' if color != "#D0D0D0" else "#666666")

    # D-odd count
    ax.text(xpos + w/2, ypos + h * 0.40, count,
            ha='center', va='center', fontsize=11,
            color='white' if color != "#D0D0D0" else "#888888")

    # Symbol
    symcolor = '#2ECC40' if sym == '✓' else '#FF4136'
    ax.text(xpos + w/2, ypos + h * 0.18, sym,
            ha='center', va='center', fontsize=20, fontweight='bold',
            color=symcolor)

# Axis labels
ax.text(x0 + w/2, y0 + 2*h + 0.06, "Y: low (8–11)\na = 0",
        ha='center', va='center', fontsize=11, color='#333')
ax.text(x0 + 3*w/2, y0 + 2*h + 0.06, "Y: high (12–15)\na = 1",
        ha='center', va='center', fontsize=11, color='#333')

ax.text(x0 - 0.07, y0 + 3*h/2, "X: low\n(8–11)\na = 0",
        ha='center', va='center', fontsize=11, color='#333', rotation=0)
ax.text(x0 - 0.07, y0 + h/2, "X: high\n(12–15)\na = 1",
        ha='center', va='center', fontsize=11, color='#333', rotation=0)

# Title
ax.text(x0 + w, y0 + 2*h + 0.16,
        "Four Sectors of K = 4 Binary Multiplication",
        ha='center', va='center', fontsize=15, fontweight='bold', color='#222')

# Theorem 4 annotation
ax.annotate(
    "Theorem 4: always D-even\n$XY \\geq 9 \\cdot 2^{2K-4} > 2^{2K-1}$",
    xy=(x0 + 3*w/2, y0 + h/2),
    xytext=(x0 + 2*w + 0.05, y0 + h * 0.3),
    fontsize=9, color='#AA3333', fontstyle='italic',
    arrowprops=dict(arrowstyle='->', color='#AA3333', lw=1.5),
    ha='left', va='center',
    bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF5F5', edgecolor='#DD9999')
)

# R(K) annotation
ax.annotate(
    "$R(K) = \\sigma_{10} \\,/\\, \\sigma_{00} \\;\\to\\; -\\pi$",
    xy=(x0 + w/2, y0 + h/2),
    xytext=(x0 - 0.05, y0 - 0.08),
    fontsize=11, fontweight='bold', color='#2255AA',
    arrowprops=dict(arrowstyle='->', color='#2255AA', lw=1.5,
                    connectionstyle="arc3,rad=0.2"),
    ha='center', va='center',
    bbox=dict(boxstyle='round,pad=0.3', facecolor='#F0F5FF', edgecolor='#99AADD')
)
# Second arrow from (1,0)
ax.annotate("",
    xy=(x0 + 3*w/2, y0 + h/2),
    xytext=(x0 + 0.05, y0 - 0.04),
    arrowprops=dict(arrowstyle='->', color='#CC8800', lw=1.5,
                    connectionstyle="arc3,rad=-0.3"))

ax.set_xlim(-0.15, 1.15)
ax.set_ylim(-0.15, 1.15)
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()
plt.savefig("papers/carry-arithmetic-P1-pi-spectral/figures/fig_sector_layout.png",
            dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print("OK: fig_sector_layout.png")
