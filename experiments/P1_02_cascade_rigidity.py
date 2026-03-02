"""
P1-RP: Prove cascade value rigidity (val = ±1) in 7/9 sub-cells
=================================================================

From P1-B4: When decomposing sector (0,0) by (c₂,c₃,c₄) =
(carries[D-2], carries[D-3], carries[D-4]):

  val = -1 exactly for ALL pairs in cells: (0,1,0), (1,0,0), (1,0,1)
  val = +1 exactly for ALL pairs in cells: (0,1,2), (1,2,2)
  val =  0 (exact cancellation) in cells: (0,1,1), (1,1,*)

Only (0,0,0) and (0,0,1) have variable val.

Question: WHY is val = carries[M-1] - 1 constant in these cells?

Strategy:
  M = max position with carries[m] > 0. For D-odd pairs (product has D bits),
  the highest carry positions are near D. If c₂ = carries[D-2] ≥ 1 and
  c₃ = carries[D-3] is known, then M is determined, and carries[M-1] follows
  from the propagation structure.

  This script:
  A) Verifies the rigidity pattern for K=7..13
  B) For each rigid cell, determines WHICH M values occur
  C) Shows that M is determined by (c₂, c₃, c₄) in the rigid cells
  D) Proves carries[M-1] follows deterministically
"""

import sys
import time
from collections import defaultdict

def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def analyze_rigidity(K):
    """For each pair, compute val and the full carry profile near D."""
    D = 2*K - 1
    lo = 1 << (K-1)
    hi = 1 << K

    # Accumulators: cell_data[(c2,c3,c4)] = list of (val, M, carries[M-1], carries near D)
    cell_info = defaultdict(lambda: {'vals': defaultdict(int), 'Ms': defaultdict(int),
                                      'cM1': defaultdict(int), 'total': 0})

    for p in range(lo, hi):
        a1 = (p >> (K-2)) & 1
        for q in range(lo, hi):
            prod = p * q
            if prod.bit_length() != D:
                continue
            c1 = (q >> (K-2)) & 1

            sec = a1 * 2 + c1
            if sec != 0:
                continue

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
            cd2 = carries[D-2] if D >= 3 else 0
            cd3 = carries[D-3] if D >= 4 else 0
            cd4 = carries[D-4] if D >= 5 else 0

            key = (cd2, cd3, cd4)
            cell_info[key]['vals'][val] += 1
            cell_info[key]['Ms'][M] += 1
            cell_info[key]['cM1'][carries[M-1]] += 1
            cell_info[key]['total'] += 1

    return cell_info


flush("=" * 76)
flush("  P1-RP: CASCADE VALUE RIGIDITY ANALYSIS")
flush("=" * 76)

for K in range(7, 14):
    D = 2*K - 1
    t0 = time.time()
    flush(f"\n  K={K}, D={D}")

    info = analyze_rigidity(K)
    elapsed = time.time() - t0

    flush(f"  ({elapsed:.1f}s)  Cells found: {len(info)}")
    flush(f"  {'cell (c2,c3,c4)':>20} {'n':>8} {'vals':>30} {'M values':>30} {'c[M-1]':>20} {'rigid?':>7}")

    for key in sorted(info.keys()):
        ci = info[key]
        n = ci['total']
        val_str = str(dict(ci['vals']))
        m_str = str(dict(ci['Ms']))
        cm1_str = str(dict(ci['cM1']))

        # Check rigidity: all vals the same?
        unique_vals = list(ci['vals'].keys())
        rigid = "YES" if len(unique_vals) == 1 else ("CANCEL" if sum(v*c for v,c in ci['vals'].items()) == 0 else "no")

        flush(f"  {str(key):>20} {n:8d} {val_str:>30} {m_str:>30} {cm1_str:>20} {rigid:>7}")

    # ── Analyze the WHY for rigid cells ──────────────────────────
    if K <= 10:
        flush(f"\n  WHY analysis for K={K}:")
        for key in sorted(info.keys()):
            ci = info[key]
            unique_vals = list(ci['vals'].keys())
            unique_Ms = sorted(ci['Ms'].keys())
            unique_cM1 = sorted(ci['cM1'].keys())

            if len(unique_vals) == 1:
                v = unique_vals[0]
                flush(f"    Cell {key}: val={v} always. M ∈ {unique_Ms}, c[M-1] ∈ {unique_cM1}")
                if v == -1:
                    flush(f"      → carries[M-1] = 0 for all pairs. M is the first nonzero carry from the top.")
                    flush(f"        So the carry drops from >0 at M to 0 at M-1. Claim: c₂ or c₃ constraints force this.")
                elif v == 1:
                    flush(f"      → carries[M-1] = 2 for all pairs. M is the first nonzero carry from the top.")
                    flush(f"        So the carry rises to 2 at M-1 before dropping. This happens when c₃≥2 or c₄≥2.")
            elif sum(v*c for v,c in ci['vals'].items()) == 0:
                flush(f"    Cell {key}: EXACT CANCELLATION. M ∈ {unique_Ms}")
                flush(f"      val distribution: {dict(ci['vals'])}")

flush(f"\n{'='*76}")
flush(f"  DONE")
flush(f"{'='*76}")
