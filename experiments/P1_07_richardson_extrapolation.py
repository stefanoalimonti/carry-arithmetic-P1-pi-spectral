"""
P1_07 — Richardson extrapolation of R(K) → −π.

Computes exact R(K) = S10/S00 for small K via exhaustive enumeration,
supplements with E44/E45 high-precision data for large K, and applies
multiple acceleration methods to extrapolate R(∞).

Methods:
  1. Neville–Aitken Richardson table (Ansatz 1: expansion in 1/K²)
  2. Neville–Aitken Richardson table (Ansatz 2: expansion in ρ^K, ρ = 1/2)
  3. Iterated Aitken Δ² (Shanks transformation)

References: Paper P1 §9.2 (Proposition 3); Richardson & Gaunt (1927).
"""

from fractions import Fraction
import math
import sys

PI = math.pi

# --- Exact integer counts from E44 (K ≤ 20) and E45 (K = 21) ---
# R(K) = S10 / S00; both are exact integers.
# K ≤ 13: computed inline (see enumerate_exact below).
# K ≥ 14: from E44/E45 C enumerator output.
E44_E45_COUNTS = {
    19: (-494176763, 1536648367),
    20: (-1967962747, 6149608524),
    21: (-7852453718, 24604885342),
}

# Approximate R(K) for K=14–18 from E44 output (4-decimal precision).
# Full precision requires running E44 (C program, ~hours for K≥16).
R_APPROX = {
    14: -2.6225,
    15: -2.8241,
    16: -2.9539,
    17: -3.0348,
    18: -3.0822,
}


def enumerate_exact(K):
    """Exhaustive enumeration of cascade-val sector ratio for K-bit operands.

    Returns (S00, S10) as Python ints (exact).
    """
    lo = 1 << (K - 1)
    hi = 1 << K
    D = 2 * K - 1
    bit_D1 = 1 << (D - 1)
    bit_D = 1 << D

    S00 = 0
    S10 = 0

    for p in range(lo, hi):
        a = (p >> (K - 2)) & 1
        for q in range(lo, hi):
            prod = p * q
            if prod < bit_D1 or prod >= bit_D:
                continue
            c = (q >> (K - 2)) & 1

            carries = [0] * (D + 1)
            for pos in range(D):
                conv = 0
                for i in range(K):
                    j = pos - i
                    if 0 <= j < K:
                        conv += ((p >> i) & 1) * ((q >> j) & 1)
                total = conv + carries[pos]
                carries[pos + 1] = total >> 1

            m_stop = -1
            for pos in range(D - 1, -1, -1):
                if carries[pos] > 0:
                    m_stop = pos
                    break
            if m_stop < 0:
                continue

            val = carries[m_stop - 1] - 1 if m_stop > 0 else -1

            if a == 0 and c == 0:
                S00 += val
            elif a == 1 and c == 0:
                S10 += val

    return S00, S10


def compute_R_table(K_max_compute=13):
    """Build table of exact R(K) values.

    Computes K = 4..K_max_compute by enumeration; uses E44/E45 for larger K.
    Returns dict K -> (R_float, precision_digits).
    """
    R = {}
    print(f"Computing exact R(K) for K = 4..{K_max_compute}...")
    sys.stdout.flush()

    for K in range(4, K_max_compute + 1):
        S00, S10 = enumerate_exact(K)
        r = Fraction(S10, S00)
        R[K] = (float(r), 15)
        print(f"  K={K:2d}: S00={S00:>14d}  S10={S10:>14d}  "
              f"R = {float(r):+.12f}")
        sys.stdout.flush()

    for K, r_approx in sorted(R_APPROX.items()):
        if K not in R:
            R[K] = (r_approx, 4)
            print(f"  K={K:2d}: R = {r_approx:+.12f}  [E44 approx, 4 dec]")

    for K, (s00, s10) in sorted(E44_E45_COUNTS.items()):
        r = Fraction(s10, s00)
        R[K] = (float(r), 12)
        print(f"  K={K:2d}: S00={s00:>14d}  S10={s10:>14d}  "
              f"R = {float(r):+.12f}  [E45 exact]")

    return R


def neville_table(K_vals, R_vals, h_func):
    """Neville–Aitken extrapolation table to h = 0.

    h_func(K) defines the step-size variable (e.g., 1/K² or (1/2)^K).
    Returns the triangular table; table[j][0] is the j-th level estimate.
    """
    n = len(K_vals)
    h = [h_func(k) for k in K_vals]
    table = [list(R_vals)]
    for j in range(1, n):
        col = []
        prev = table[j - 1]
        for i in range(n - j):
            val = (h[i] * prev[i + 1] - h[i + j] * prev[i]) / (h[i] - h[i + j])
            col.append(val)
        table.append(col)
    return table


def aitken_delta2(s):
    """Aitken Δ² acceleration on a sequence."""
    out = []
    for i in range(len(s) - 2):
        denom = s[i + 2] - 2 * s[i + 1] + s[i]
        if abs(denom) < 1e-30:
            out.append(s[i + 2])
        else:
            out.append(s[i] - (s[i + 1] - s[i])**2 / denom)
    return out


def iterated_aitken(s, depth=3):
    """Apply Aitken Δ² repeatedly."""
    current = list(s)
    for _ in range(depth):
        if len(current) < 3:
            break
        current = aitken_delta2(current)
    return current


def sig_digits(est):
    err = abs(est + PI)
    if err == 0:
        return float('inf')
    return -math.log10(err / PI)


def main():
    K_max = 13
    if len(sys.argv) > 1:
        K_max = int(sys.argv[1])

    R_table = compute_R_table(K_max)

    Ks = sorted(R_table.keys())
    Rs = [R_table[k][0] for k in Ks]

    print()
    print("=" * 70)
    print("Richardson Extrapolation of R(K) → −π")
    print("=" * 70)
    print()

    print(f"{'K':>3s}  {'R(K)':>16s}  {'|R(K)+π|':>12s}  {'sig.dig':>8s}  {'g(K)/g(K-1)':>12s}")
    prev_gap = None
    for K, R in zip(Ks, Rs):
        gap = R + PI
        sd = sig_digits(R)
        ratio_str = ""
        if prev_gap is not None and abs(prev_gap) > 1e-15:
            ratio_str = f"{gap / prev_gap:12.3f}"
        print(f"{K:3d}  {R:+16.12f}  {abs(gap):12.6e}  {sd:8.1f}  {ratio_str}")
        prev_gap = gap
    print()

    # Use odd K ≥ 7 for Richardson (matches paper)
    K_odd = [k for k in Ks if k >= 7 and k % 2 == 1]
    R_odd = [R_table[k][0] for k in K_odd]

    # --- Method 1: Neville with h = 1/K² (Ansatz 1) ---
    print("Method 1: Richardson–Neville, h = 1/K² (polynomial ansatz)")
    table1 = neville_table(K_odd, R_odd, lambda k: 1.0 / k**2)
    for j in range(len(table1)):
        est = table1[j][0]
        sd = sig_digits(est)
        tag = " ← best" if j == len(table1) - 1 else ""
        print(f"  Level {j:2d}: R_∞ = {est:+16.12f}  sig.dig = {sd:5.1f}{tag}")
    print()

    # --- Method 2: Neville with h = (1/2)^K (exponential ansatz) ---
    print("Method 2: Richardson–Neville, h = (1/2)^K (exponential ansatz)")
    table2 = neville_table(K_odd, R_odd, lambda k: 0.5**k)
    for j in range(len(table2)):
        est = table2[j][0]
        sd = sig_digits(est)
        tag = " ← best" if j == len(table2) - 1 else ""
        print(f"  Level {j:2d}: R_∞ = {est:+16.12f}  sig.dig = {sd:5.1f}{tag}")
    print()

    # --- Method 3: Iterated Aitken Δ² ---
    print("Method 3: Iterated Aitken Δ² (all K ≥ 7)")
    K_tail = [k for k in Ks if k >= 7]
    R_tail = [R_table[k][0] for k in K_tail]

    current = list(R_tail)
    for depth in range(1, 5):
        if len(current) < 3:
            break
        current = aitken_delta2(current)
        best = current[-1] if current else R_tail[-1]
        sd = sig_digits(best)
        print(f"  Depth {depth}: R_∞ = {best:+16.12f}  sig.dig = {sd:5.1f}"
              f"  (from {len(current)} values)")
    print()

    # --- Method 2b: Neville exponential on high-K tail ---
    if len(K_odd) >= 5:
        K_tail5 = K_odd[-5:]
        R_tail5 = [R_table[k][0] for k in K_tail5]
        print(f"Method 2b: Richardson–Neville exponential, K = {K_tail5}")
        table2b = neville_table(K_tail5, R_tail5, lambda k: 0.5**k)
        for j in range(len(table2b)):
            est = table2b[j][0]
            sd = sig_digits(est)
            tag = " ← best" if j == len(table2b) - 1 else ""
            print(f"  Level {j:2d}: R_∞ = {est:+16.12f}  sig.dig = {sd:5.1f}{tag}")
        print()

    # --- Summary ---
    candidates = []
    if table1[-1]:
        candidates.append(("Polynomial Richardson", table1[-1][0]))
    if table2[-1]:
        candidates.append(("Exponential Richardson", table2[-1][0]))
    if current:
        candidates.append(("Iterated Aitken", current[-1]))

    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Target: −π = {-PI:.12f}")
    print()
    for name, est in candidates:
        sd = sig_digits(est)
        print(f"  {name:30s}: R_∞ = {est:+.12f}  sig.digits = {sd:.1f}")

    best_est = max(candidates, key=lambda x: sig_digits(x[1]))
    print()
    print(f"  Best: {best_est[0]} → {sig_digits(best_est[1]):.1f} significant digits")

    # --- LOO cross-validation for best method ---
    print()
    print("Leave-one-out (exponential Richardson):")
    for idx in range(len(K_odd)):
        K_sub = K_odd[:idx] + K_odd[idx + 1:]
        R_sub = [R_table[k][0] for k in K_sub]
        t = neville_table(K_sub, R_sub, lambda k: 0.5**k)
        est = t[-1][0]
        sd = sig_digits(est)
        print(f"  Leave out K={K_odd[idx]:2d}: R_∞ = {est:+.12f}  "
              f"sig.digits = {sd:.1f}")


if __name__ == "__main__":
    main()
