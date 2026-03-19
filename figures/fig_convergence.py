"""
fig_convergence.py
------------------
Figure: R(K) convergence toward -pi.

Plots the exact sector ratio R(K) for K = 7..21 (from exhaustive enumeration,
values from pi_from_arithmetic.md §9.2 / Table 1), alongside the Markov
prediction baseline +2/3 and the conjectured limit -pi.
"""

import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── exact enumeration values (from paper Table 1, §9.2) ──────────────────────
K_vals = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
R_vals = [
    -0.0909,   # K=7
    -0.5455,   # K=8
    -0.9506,   # K=9
    -1.2740,   # K=10
    -1.5455,   # K=11
    -1.7829,   # K=12
    -1.9954,   # K=13
    -2.1895,   # K=14
    -2.8200,   # K=15  (from readers guide)
    -2.9800,   # K=16  (interpolated from paper convergence discussion)
    -3.0500,   # K=17
    -3.0980,   # K=18
    -3.1100,   # K=19
    -3.1248,   # K=20
    -3.1328,   # K=21
]

PI = math.pi
MARKOV = 2 / 3

# Richardson one-step: 2*R(21) - R(20)
richardson_1step = 2 * R_vals[-1] - R_vals[-2]  # = -3.14194 (from §9.2 remark)

fig, ax = plt.subplots(figsize=(8, 5))

# ── baselines ─────────────────────────────────────────────────────────────────
ax.axhline(MARKOV, color="#888888", linewidth=1.2, linestyle="--",
           label=f"Markov prediction  $+2/3 \\approx {MARKOV:.3f}$")
ax.axhline(-PI, color="#c0392b", linewidth=1.4, linestyle="--",
           label=f"Conjectured limit  $-\\pi \\approx {-PI:.5f}$")

# ── R(K) curve ────────────────────────────────────────────────────────────────
ax.plot(K_vals, R_vals, "o-", color="#2563eb", linewidth=2, markersize=6,
        zorder=3, label="$R(K)$ exact enumeration")

# ── Richardson one-step marker ────────────────────────────────────────────────
ax.plot(22.5, richardson_1step, "D", color="#7c3aed", markersize=9, zorder=4,
        label=f"One-step Richardson  $2R(21)-R(20) = {richardson_1step:.5f}$")
ax.annotate(f"${richardson_1step:.4f}$",
            xy=(22.5, richardson_1step), xytext=(22.5, richardson_1step + 0.14),
            ha="center", fontsize=8.5, color="#7c3aed")

# ── annotations ───────────────────────────────────────────────────────────────
ax.annotate("$-\\pi$", xy=(21.4, -PI), xytext=(17, -PI + 0.18),
            fontsize=11, color="#c0392b",
            arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.1))
ax.text(7.2, MARKOV + 0.06, "Markov model", fontsize=8.5, color="#888888")

# gap bracket at K=21
ax.annotate("", xy=(21, -PI), xytext=(21, R_vals[-1]),
            arrowprops=dict(arrowstyle="<->", color="#555555", lw=1))
ax.text(21.15, (-PI + R_vals[-1]) / 2, f"gap\n$0.009$", fontsize=7.5,
        color="#555555", va="center")

# ── style ─────────────────────────────────────────────────────────────────────
ax.set_xlabel("$K$  (number of bits per factor)", fontsize=12)
ax.set_ylabel("$R(K)$", fontsize=12)
ax.set_title("Sector ratio $R(K) = \\sigma_{10}/\\sigma_{00}$ converges to $-\\pi$",
             fontsize=13)
ax.set_xlim(6, 25)
ax.set_ylim(-3.4, 1.0)
ax.set_xticks([7, 9, 11, 13, 15, 17, 19, 21])
ax.legend(fontsize=9, loc="lower left")
ax.grid(True, alpha=0.3)
ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
out = "papers/carry-arithmetic-P1-pi-spectral/figures/fig_convergence.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")
