"""
P1-B4: Test the Bernoulli hierarchy Level 4 prediction |B₄| = 1/30
===================================================================

Key insight from run 1: The D-3 level encodes B₃=0 (σ(c2=1,c3=1) ≡ 0).
To find |B₄|=1/30, we likely need the carry at D-4.

This version decomposes by c₂=carries[D-2], c₃=carries[D-3], c₄=carries[D-4]
and tests all plausible ratio candidates against 1/30.

The Euler-Maclaurin hierarchy suggests nested sub-partitions:
  Level 2: α₂ = σ^{c₂≥1}/σ^{c₂=0}  → B₂ = 1/6
  Level 3: σ(c₂=1,c₃=1) ≡ 0          → B₃ = 0  (confirmed!)
  Level 4: α₄ = ??? within c₂=1,c₃=0 sub-population → |B₄| = 1/30 ?
"""

import sys
import time
from collections import defaultdict
from math import pi

def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def compute_data(K):
    """Exact enumeration for K-bit inputs."""
    D = 2*K - 1
    lo = 1 << (K-1)
    hi = 1 << K

    S00_c0 = 0
    S00_c1 = 0
    n00_c0 = 0
    n00_c1 = 0

    # 3D accumulators: by (c_{D-2}, c_{D-3}, c_{D-4})
    sigma5 = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    count5 = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    for p in range(lo, hi):
        a1 = (p >> (K-2)) & 1
        for q in range(lo, hi):
            prod = p * q
            if prod.bit_length() != D:
                continue
            c1 = (q >> (K-2)) & 1

            carries = [0] * (D + 1)
            for d in range(D):
                cv = 0
                for i in range(max(0, d - K + 1), min(d, K - 1) + 1):
                    cv += ((p >> i) & 1) * ((q >> (d - i)) & 1)
                carries[d + 1] = (cv + carries[d]) >> 1

            M = None
            for m in range(D, 0, -1):
                if carries[m] > 0:
                    M = m
                    break
            if M is None:
                continue

            val = carries[M - 1] - 1

            sec = a1 * 2 + c1
            if sec != 0:
                continue

            cd2 = carries[D-2] if D >= 3 else 0
            cd3 = carries[D-3] if D >= 4 else 0
            cd4 = carries[D-4] if D >= 5 else 0

            if cd2 == 0:
                S00_c0 += val
                n00_c0 += 1
            else:
                S00_c1 += val
                n00_c1 += 1

            sigma5[cd2][cd3][cd4] += val
            count5[cd2][cd3][cd4] += 1

    return {
        'S00_c0': S00_c0, 'S00_c1': S00_c1,
        'n00_c0': n00_c0, 'n00_c1': n00_c1,
        'sigma5': sigma5, 'count5': count5, 'D': D,
    }


flush("=" * 76)
flush("  P1-B4: BERNOULLI HIERARCHY — LEVELS 2, 3, 4")
flush("=" * 76)
flush(f"\n  Level 2 target: α₂ → 1/6 = {1/6:.10f}")
flush(f"  Level 3 signal: B₃ = 0 → σ(c₂=1,c₃=1) ≡ 0")
flush(f"  Level 4 target: α₄ → 1/30 = {1/30:.10f}")

results = []

for K in range(7, 22):
    D = 2*K - 1
    lo = 1 << (K-1)
    hi = 1 << K
    n_pairs = (hi - lo) ** 2

    t0 = time.time()
    flush(f"\n  K={K:2d}  D={D:2d}  pairs={n_pairs:>12,d} ... ", end="")

    data = compute_data(K)
    elapsed = time.time() - t0

    alpha = data['S00_c1'] / data['S00_c0'] if data['S00_c0'] != 0 else float('nan')

    s5 = data['sigma5']
    c5 = data['count5']

    flush(f"({elapsed:6.1f}s)  α₂={alpha:+.8f}")

    # Print full 3D distribution
    flush(f"    c₂ \\ c₃ \\ c₄ distribution (σ values):")
    for cd2 in sorted(s5.keys()):
        for cd3 in sorted(s5[cd2].keys()):
            cd4_vals = sorted(s5[cd2][cd3].keys())
            parts = ", ".join(f"c₄={cd4}: σ={s5[cd2][cd3][cd4]:>9d} n={c5[cd2][cd3][cd4]:>9d}"
                              for cd4 in cd4_vals)
            flush(f"      c₂={cd2}, c₃={cd3}: {parts}")

    # ── Candidate ratios for Level 4 ─────────────────────────────
    # The prediction says "sub-partition of the c=1 population"
    # c=1 means c_{D-2}=1. Within that: c₃ ∈ {0, 2} (c₃=1 is zero).
    # Further decompose by c₄:

    # (a) Within c₂=1, c₃=0: ratio by c₄
    #    α₄ᵃ = σ(1,0,c₄≥1) / σ(1,0,0)
    s_1_0_0 = s5.get(1, {}).get(0, {}).get(0, 0)
    s_1_0_ge1 = sum(s5.get(1, {}).get(0, {}).get(cd4, 0)
                    for cd4 in s5.get(1, {}).get(0, {}) if cd4 >= 1)
    ratio_a = s_1_0_ge1 / s_1_0_0 if s_1_0_0 != 0 else float('nan')

    # (b) Within c₂=1, c₃=2: ratio by c₄
    s_1_2_0 = s5.get(1, {}).get(2, {}).get(0, 0)
    s_1_2_ge1 = sum(s5.get(1, {}).get(2, {}).get(cd4, 0)
                    for cd4 in s5.get(1, {}).get(2, {}) if cd4 >= 1)
    ratio_b = s_1_2_ge1 / s_1_2_0 if s_1_2_0 != 0 else float('nan')

    # (c) Within c₂=0, c₃=0: ratio by c₄
    s_0_0_0 = s5.get(0, {}).get(0, {}).get(0, 0)
    s_0_0_ge1 = sum(s5.get(0, {}).get(0, {}).get(cd4, 0)
                    for cd4 in s5.get(0, {}).get(0, {}) if cd4 >= 1)
    ratio_c = s_0_0_ge1 / s_0_0_0 if s_0_0_0 != 0 else float('nan')

    # (d) Within c₂=0, c₃=1: ratio by c₄
    s_0_1_0 = s5.get(0, {}).get(1, {}).get(0, 0)
    s_0_1_ge1 = sum(s5.get(0, {}).get(1, {}).get(cd4, 0)
                    for cd4 in s5.get(0, {}).get(1, {}) if cd4 >= 1)
    ratio_d = s_0_1_ge1 / s_0_1_0 if s_0_1_0 != 0 else float('nan')

    # (e) Overall c₄ ratio within c₂=1:
    #    σ^{c₂=1,c₄≥1} / σ^{c₂=1,c₄=0}
    s_c2eq1_c4_0 = sum(s5.get(1, {}).get(cd3, {}).get(0, 0)
                       for cd3 in s5.get(1, {}))
    s_c2eq1_c4_ge1 = sum(s5.get(1, {}).get(cd3, {}).get(cd4, 0)
                         for cd3 in s5.get(1, {})
                         for cd4 in s5.get(1, {}).get(cd3, {}) if cd4 >= 1)
    ratio_e = s_c2eq1_c4_ge1 / s_c2eq1_c4_0 if s_c2eq1_c4_0 != 0 else float('nan')

    # (f) σ(1,0,1) / σ(1,0,0) — the most specific
    s_1_0_1 = s5.get(1, {}).get(0, {}).get(1, 0)
    ratio_f = s_1_0_1 / s_1_0_0 if s_1_0_0 != 0 else float('nan')

    # (g) σ(0,0,1) / σ(0,0,0) — nested α within c₂=0,c₃=0
    s_0_0_1 = s5.get(0, {}).get(0, {}).get(1, 0)
    ratio_g = s_0_0_1 / s_0_0_0 if s_0_0_0 != 0 else float('nan')

    results.append({
        'K': K, 'D': D, 'alpha': alpha, 'time': elapsed,
        'a': ratio_a, 'b': ratio_b, 'c': ratio_c, 'd': ratio_d,
        'e': ratio_e, 'f': ratio_f, 'g': ratio_g,
    })

    flush(f"    Ratios: a(1,0,≥1)/(1,0,0)={ratio_a:+.6f}  "
          f"b(1,2,≥1)/(1,2,0)={ratio_b:+.6f}  "
          f"f(1,0,1)/(1,0,0)={ratio_f:+.6f}  "
          f"g(0,0,1)/(0,0,0)={ratio_g:+.6f}")

    if elapsed > 180:
        flush(f"  (stopping — next K would be too slow)")
        break


# ── SUMMARY ───────────────────────────────────────────────────────
flush(f"\n{'=' * 76}")
flush(f"  CONVERGENCE SUMMARY")
flush(f"{'=' * 76}")

target = 1/30

flush(f"\n  Level 2: α₂ → 1/6 = {1/6:.10f}")
flush(f"  {'K':>3s}  {'α₂':>12s}  {'|α₂-1/6|':>12s}")
for r in results:
    gap = abs(r['alpha'] - 1/6)
    flush(f"  {r['K']:3d}  {r['alpha']:12.8f}  {gap:12.8f}")

for label, key in [
    ('a: σ(1,0,c₄≥1)/σ(1,0,0) — within c₂=1,c₃=0 by c₄', 'a'),
    ('b: σ(1,2,c₄≥1)/σ(1,2,0) — within c₂=1,c₃=2 by c₄', 'b'),
    ('c: σ(0,0,c₄≥1)/σ(0,0,0) — within c₂=0,c₃=0 by c₄', 'c'),
    ('d: σ(0,1,c₄≥1)/σ(0,1,0) — within c₂=0,c₃=1 by c₄', 'd'),
    ('e: σ(c₂=1,c₄≥1)/σ(c₂=1,c₄=0) — overall in c₂=1', 'e'),
    ('f: σ(1,0,1)/σ(1,0,0) — specific cell', 'f'),
    ('g: σ(0,0,1)/σ(0,0,0) — nested α', 'g'),
]:
    flush(f"\n  Candidate {label}")
    flush(f"  → 1/30 = {target:.10f}?")
    flush(f"  {'K':>3s}  {'ratio':>12s}  {'|r-1/30|':>12s}  {'trend':>8s}")
    prev_gap = None
    for r in results:
        v = r[key]
        if v != v or abs(v) > 1e10:
            flush(f"  {r['K']:3d}  {'NaN/Inf':>12s}")
            continue
        gap = abs(v - target)
        trend = ""
        if prev_gap is not None:
            trend = "↓" if gap < prev_gap else "↑"
        flush(f"  {r['K']:3d}  {v:12.8f}  {gap:12.8f}  {trend:>8s}")
        prev_gap = gap

flush(f"\n{'=' * 76}")
flush(f"  DONE")
flush(f"{'=' * 76}")
