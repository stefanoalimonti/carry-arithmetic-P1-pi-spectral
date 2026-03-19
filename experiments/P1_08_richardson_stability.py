"""
P1_08 — One-step Richardson check and topological stability of minima.

Two independent verifications:

Part A — One-step Richardson:
  Applies the single deferred-limit step  2*R(21) - R(20)  to the exact
  sector-ratio data, yielding a 3.5-significant-digit estimate of R(inf).
  This is a stripped-down consistency check of the full Neville-Aitken
  result of P1_07: it confirms that 3+ digit precision does not require
  the full extrapolation table.

Part B — Positional stability of |mu_chi4| minima:
  Compares the K = 21 and K = 999 profiles of |mu_chi4(1/2 + it)| on a
  grid with spacing dt = 0.05.  For each of 15 minima near the known
  zeros of L(s, chi4), records the positional shift.  The expected result
  is shift = 0.000 on all minima (see also L22 in paper L).

References: Paper P1 §9.2 (Proposition 3, remark), §9.5.
"""

import math
import sys
import os
import numpy as np

PI = math.pi

# ── Exact sector-ratio data (E44/E45 enumeration) ─────────────────────────
# Source: E44 for K <= 20, E45 for K = 21.
R_EXACT = {
    4:  +0.2500,
    5:  +0.2500,
    6:  +0.09302,
    7:  -0.09149,
    8:  -0.38480,
    9:  -0.72617,
    10: -1.11837,
    11: -1.55346,
    12: -1.96958,
    13: -2.33876,
    14: -2.62253,
    15: -2.82406,
    16: -2.95386,
    17: -3.03482,
    18: -3.08216,
    19: -3.10951,
    20: -3.12486,
    21: -3.13340,
}


def main():
    print("=" * 70)
    print("P1_08: ONE-STEP RICHARDSON AND STABILITY OF MINIMA POSITIONS")
    print("=" * 70)

    # ── Part A: One-step Richardson ───────────────────────────────────────
    print("\n--- Part A: One-step Richardson  2*R(21) - R(20) ---")
    R20 = R_EXACT[20]
    R21 = R_EXACT[21]
    R_rich = 2 * R21 - R20
    gap = abs(R_rich + PI)
    sig_digits = -math.log10(gap / PI) if gap > 0 else float('inf')

    print(f"\n  R(20) = {R20:.8f}")
    print(f"  R(21) = {R21:.8f}")
    print(f"  2*R(21) - R(20) = {R_rich:.8f}")
    print(f"  -pi              = {-PI:.8f}")
    print(f"  gap              = {gap:.2e}  ({sig_digits:.1f} significant digits)")

    # Also compute gap ratios for K = 19..21 to confirm asymptotic regime
    print("\n  Gap ratios g(K)/g(K-1) for recent K:")
    print(f"  {'K':>4}  {'R(K)':>12}  {'g(K)=|R+pi|':>14}  {'g/g_prev':>10}")
    prev_g = None
    for K in [18, 19, 20, 21]:
        g = abs(R_EXACT[K] + PI)
        ratio = g / prev_g if prev_g else float('nan')
        print(f"  {K:>4d}  {R_EXACT[K]:>12.6f}  {g:>14.6f}  {ratio:>10.4f}")
        prev_g = g

    print(f"\n  Convergence rate -> 1/2 (Diaconis-Fulman spectral radius)")

    # ── Part B: Positional stability of |mu_chi4| minima ─────────────────
    print("\n--- Part B: Positional stability of |mu_chi4(1/2+it)| minima ---")

    # Import carry-bank library from paper L
    paper_L_exp = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "carry-arithmetic-L-dirichlet-bridge", "experiments"
    )
    sys.path.insert(0, os.path.abspath(paper_L_exp))

    try:
        from _shared import build_highk_bank, mu_chi_from_record, primitive_characters
    except ImportError as e:
        print(f"  (skipping Part B: cannot import carry bank — {e})")
        return

    bank = build_highk_bank()
    chi4_fn = {c["name"]: c for c in primitive_characters()}["chi4"]["fn"]

    ZEROS_L = [6.0209, 10.2438, 12.9881, 16.3426, 18.2920,
               21.4506, 23.2784, 25.7288, 28.3596, 29.6564,
               32.5922, 34.2000, 36.1429, 38.5119, 40.3227]

    sigma = 0.5
    dt = 0.05
    ts = np.arange(1.0, 55.0, dt)

    print(f"\n  Computing profiles for K = 21 and K = 999 (grid dt = {dt}) ...")
    mu21  = np.array([mu_chi_from_record(bank[21],  chi4_fn, complex(sigma, t)) for t in ts])
    mu999 = np.array([mu_chi_from_record(bank[999], chi4_fn, complex(sigma, t)) for t in ts])

    def find_min(arr, t_c, win=1.5):
        mask = (ts >= t_c - win) & (ts <= t_c + win)
        if not np.any(mask): return float('nan')
        return ts[mask][np.argmin(np.abs(arr[mask]))]

    print(f"\n  {'#':>3}  {'t_zero_L':>10}  {'t_min K=21':>12}  {'t_min K=999':>12}  {'shift':>8}")
    shifts = []
    for i, tz in enumerate(ZEROS_L):
        t21  = find_min(mu21,  tz)
        t999 = find_min(mu999, tz)
        shift = abs(t21 - t999)
        shifts.append(shift)
        print(f"  {i+1:>3d}  {tz:>10.4f}  {t21:>12.4f}  {t999:>12.4f}  {shift:>8.4f}")

    print(f"\n  Mean shift: {np.mean(shifts):.4f}  ({np.mean(shifts)/dt:.1f} grid points)")
    print(f"  Max  shift: {np.max(shifts):.4f}  ({np.max(shifts)/dt:.1f} grid points)")

    print("\n--- Summary ---")
    print(f"  One-step Richardson:  2R(21)-R(20) = {R_rich:.5f}  "
          f"(gap {gap:.1e} from -pi, {sig_digits:.1f} sig. digits)")
    if all(s < dt for s in shifts):
        print(f"  Positional stability: all {len(ZEROS_L)} minima have shift = 0 grid points.")
    else:
        print(f"  Positional stability: max shift = {max(shifts):.4f}")


if __name__ == "__main__":
    main()
