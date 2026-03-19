"""
fig_markov_gap.py
-----------------
Figure: The Markov prediction vs observed sector ratio.

A side-by-side comparison showing four ways the Markov model fails:
  sign, magnitude, algebraic nature, and value.
Also shows the convergence trajectory as a small inset.
"""

import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

PI = math.pi
MARKOV = 2 / 3
OBSERVED = -PI

fig, axes = plt.subplots(1, 2, figsize=(11, 5),
                         gridspec_kw={"width_ratios": [1.2, 1]})

# ── LEFT: bar comparison ──────────────────────────────────────────────────────
ax = axes[0]

labels = ["Markov\nprediction", "Observed\nlimit"]
values = [MARKOV, OBSERVED]
colors = ["#6b7280", "#c0392b"]

bars = ax.bar(labels, values, color=colors, width=0.45, zorder=3,
              edgecolor="white", linewidth=1.5)

# value annotations inside bars
ax.text(0, MARKOV / 2, f"$+2/3$\n$\\approx {MARKOV:.4f}$",
        ha="center", va="center", fontsize=12, color="white", fontweight="bold")
ax.text(1, OBSERVED / 2, f"$-\\pi$\n$\\approx {OBSERVED:.4f}$",
        ha="center", va="center", fontsize=12, color="white", fontweight="bold")

# four-way failure annotations
failure_props = dict(boxstyle="round,pad=0.3", facecolor="#fef9c3",
                     edgecolor="#d97706", alpha=0.95)
ax.text(0.5, 0.55,
        "① Sign: $+$ vs $-$\n"
        "② Magnitude: $0.67$ vs $3.14$\n"
        "③ Rational vs transcendental\n"
        "④ Finite vs irrational",
        transform=ax.transAxes, fontsize=9.5,
        va="top", ha="center", bbox=failure_props)

ax.axhline(0, color="black", lw=0.8)
ax.set_ylim(-3.8, 1.3)
ax.set_ylabel("$R(\\infty)$", fontsize=12)
ax.set_title("The Markov model fails in every way", fontsize=12)
ax.grid(axis="y", alpha=0.3)
ax.spines[["top", "right"]].set_visible(False)

# ── RIGHT: convergence trajectory ─────────────────────────────────────────────
ax2 = axes[1]

K_vals = [7, 9, 11, 13, 15, 17, 19, 20, 21]
R_vals = [-0.0909, -0.9506, -1.5455, -1.9954, -2.8200, -3.0500, -3.1100, -3.1248, -3.1328]

ax2.axhline(MARKOV, color="#6b7280", lw=1.2, ls="--", label=f"Markov $+2/3$")
ax2.axhline(-PI, color="#c0392b", lw=1.4, ls="--", label=f"Limit $-\\pi$")
ax2.plot(K_vals, R_vals, "o-", color="#2563eb", lw=2, ms=6, zorder=3,
         label="$R(K)$ exact")

# Richardson marker
ax2.plot(22.5, -3.14194, "D", color="#7c3aed", ms=9, zorder=4,
         label="Richardson $= -3.14194$")

# shade the gap
ax2.fill_between([6, 24], [MARKOV, MARKOV], [-PI, -PI],
                 alpha=0.04, color="#c0392b")

ax2.set_xlabel("$K$ (bits)", fontsize=11)
ax2.set_ylabel("$R(K)$", fontsize=11)
ax2.set_title("Convergence to $-\\pi$\n(1 trillion pairs at $K=21$)", fontsize=11)
ax2.set_xlim(6, 25)
ax2.set_ylim(-3.6, 1.0)
ax2.set_xticks([7, 9, 11, 13, 15, 17, 19, 21])
ax2.legend(fontsize=8, loc="lower left")
ax2.grid(alpha=0.3)
ax2.spines[["top", "right"]].set_visible(False)

plt.suptitle(
    "Binary multiplication carries: from rational prediction to transcendental reality",
    fontsize=13, y=1.01)
plt.tight_layout()
out = "papers/carry-arithmetic-P1-pi-spectral/figures/fig_markov_gap.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")
