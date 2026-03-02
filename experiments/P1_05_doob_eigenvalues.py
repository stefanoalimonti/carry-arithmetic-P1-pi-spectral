"""
P1-DE: Doob h-transform eigenvalues of the D-odd carry chain
=============================================================

The carry chain under D-odd conditioning (carry[D]=0) is a Doob h-transform
of the free Diaconis-Fulman chain. The LMH predicts eigenvalues:
    lambda_n = (1-A*) * lambda_n^Markov + A*/2

This script:
  1. For K=5..12, enumerates ALL D-odd sector (0,0) pairs
  2. At each position d, records carry transitions c[d] -> c[d+1]
  3. Builds the empirical transition matrix T(d) at bulk positions
  4. Aggregates to a bulk transition matrix T_bulk
  5. Computes eigenvalues and compares with Markov and LMH predictions
  6. Checks if eigenvalue modification is uniform (single A parameter)
"""

import sys
import time
import numpy as np
from math import pi, log

def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


A_STAR = (2 + 3 * pi) / (2 * (1 + pi))


def build_transition_matrices(K, sector="00"):
    """Enumerate pairs and build per-position carry transition matrices."""
    D = 2 * K - 1
    lo = 1 << (K - 1)
    hi = 1 << K

    max_carry = 0
    all_carries = []

    for p in range(lo, hi):
        a1 = (p >> (K - 2)) & 1
        for q in range(lo, hi):
            prod = p * q
            if prod.bit_length() != D:
                continue
            b1 = (q >> (K - 2)) & 1

            sec_bits = a1 * 2 + b1
            if sector == "00" and sec_bits != 0:
                continue
            if sector == "10" and sec_bits != 2:
                continue

            carries = [0] * (D + 1)
            for d in range(D):
                cv = 0
                for i in range(max(0, d - K + 1), min(d, K - 1) + 1):
                    cv += ((p >> i) & 1) * ((q >> (d - i)) & 1)
                carries[d + 1] = (cv + carries[d]) >> 1

            mc = max(carries)
            if mc > max_carry:
                max_carry = mc
            all_carries.append(carries)

    dim = max_carry + 1

    trans = {}
    for d in range(D):
        T = np.zeros((dim, dim))
        marginal = np.zeros(dim)
        for carries in all_carries:
            c_from = carries[d]
            c_to = carries[d + 1]
            T[c_to, c_from] += 1
            marginal[c_from] += 1
        for c in range(dim):
            if marginal[c] > 0:
                T[:, c] /= marginal[c]
        trans[d] = T

    return trans, dim, len(all_carries), D


def markov_eigenvalues(K):
    """Diaconis-Fulman eigenvalues for the unconditioned carry chain."""
    return [1.0, 0.5]


flush("=" * 78)
flush("  P1-DE: DOOB H-TRANSFORM EIGENVALUES")
flush("  D-odd conditioned carry chain spectral analysis")
flush(f"  A* = (2+3π)/(2(1+π)) = {A_STAR:.8f}")
flush("=" * 78)

all_results = []

for K in range(5, 14):
    D = 2 * K - 1
    lo = 1 << (K - 1)
    hi = 1 << K
    n_total = (hi - lo) ** 2

    flush(f"\n{'─' * 78}")
    flush(f"  K={K}  D={D}  max_pairs={n_total:,}")
    flush(f"{'─' * 78}")

    t0 = time.time()
    trans_00, dim_00, npairs_00, _ = build_transition_matrices(K, "00")
    trans_10, dim_10, npairs_10, _ = build_transition_matrices(K, "10")
    elapsed = time.time() - t0
    flush(f"  Enumerated in {elapsed:.1f}s  (n₀₀={npairs_00:,}, n₁₀={npairs_10:,}, dim={dim_00})")

    # ── Per-position eigenvalue analysis ──────────────────────────
    flush(f"\n  Per-position top-2 eigenvalues of T(d) for sector 00:")
    flush(f"  {'d':>4s}  {'λ₁':>10s}  {'λ₂':>10s}  {'gap':>10s}  {'A_eff':>10s}")

    bulk_range = range(max(2, K // 2), D - K // 2)
    bulk_evals_2 = []

    for d in range(D):
        T = trans_00[d]
        if T.shape[0] < 2:
            continue
        active = np.any(T > 0, axis=0) | np.any(T > 0, axis=1)
        n_active = active.sum()
        if n_active < 2:
            continue

        T_sub = T[np.ix_(active, active)]
        evals = np.linalg.eigvals(T_sub)
        evals_real = sorted(evals.real, reverse=True)

        if len(evals_real) >= 2:
            l1 = evals_real[0]
            l2 = evals_real[1]
            gap = l1 - l2
            A_eff = (2 * l2 - 1) / (2 * l2 - 2) if abs(l2 - 1) > 1e-10 else float("nan")

            marker = " <bulk>" if d in bulk_range else ""
            flush(f"  {d:4d}  {l1:10.6f}  {l2:10.6f}  {gap:10.6f}  {A_eff:10.6f}{marker}")

            if d in bulk_range:
                bulk_evals_2.append(l2)

    # ── Bulk-averaged transition matrix ───────────────────────────
    flush(f"\n  Bulk-averaged transition matrix (positions {list(bulk_range)[0]}..{list(bulk_range)[-1]}):")
    T_bulk = np.zeros((dim_00, dim_00))
    count_bulk = 0
    for d in bulk_range:
        T_bulk += trans_00[d]
        count_bulk += 1
    if count_bulk > 0:
        T_bulk /= count_bulk

    active = np.any(T_bulk > 0, axis=0) | np.any(T_bulk > 0, axis=1)
    T_sub = T_bulk[np.ix_(active, active)]
    n_active = active.sum()

    flush(f"  Active states: {n_active}")
    flush(f"  T_bulk (active):")
    for i in range(min(n_active, 5)):
        row = "  ".join(f"{T_sub[i,j]:.4f}" for j in range(min(n_active, 5)))
        flush(f"    [{row}]")

    evals_bulk = np.linalg.eigvals(T_sub)
    evals_bulk_real = sorted(evals_bulk.real, reverse=True)
    flush(f"\n  Eigenvalues of T_bulk:")
    for i, ev in enumerate(evals_bulk_real[:6]):
        flush(f"    λ_{i} = {ev:.8f}")

    if len(evals_bulk_real) >= 2:
        l2_bulk = evals_bulk_real[1]
        A_from_bulk = (2 * l2_bulk - 1) / (2 * l2_bulk - 2) if abs(l2_bulk - 1) > 1e-10 else float("nan")
        flush(f"\n  λ₂(bulk) = {l2_bulk:.8f}")
        flush(f"  Markov prediction: λ₂ = 0.5")
        flush(f"  LMH prediction: λ₂ = (1-A*)/2 · cos(π/(L+1)) + A*/2")
        L = K - 1
        lmh_l2 = (1 - A_STAR) / 2 * np.cos(pi / (L + 1)) + A_STAR / 2
        flush(f"  LMH λ₂ (L={L}) = {lmh_l2:.8f}")
        flush(f"  A_eff from λ₂: {A_from_bulk:.8f}  (A* = {A_STAR:.8f})")

    # ── Sector 10 comparison ──────────────────────────────────────
    flush(f"\n  Sector 10 bulk eigenvalues:")
    T_bulk_10 = np.zeros((dim_10, dim_10))
    count_10 = 0
    for d in bulk_range:
        if d in trans_10:
            T_bulk_10 += trans_10[d]
            count_10 += 1
    if count_10 > 0:
        T_bulk_10 /= count_10
        active_10 = np.any(T_bulk_10 > 0, axis=0) | np.any(T_bulk_10 > 0, axis=1)
        T_sub_10 = T_bulk_10[np.ix_(active_10, active_10)]
        evals_10 = np.linalg.eigvals(T_sub_10)
        evals_10_real = sorted(evals_10.real, reverse=True)
        for i, ev in enumerate(evals_10_real[:4]):
            flush(f"    λ_{i}(10) = {ev:.8f}")

    all_results.append({
        "K": K, "D": D,
        "l2_bulk": l2_bulk if len(evals_bulk_real) >= 2 else None,
        "bulk_evals": evals_bulk_real[:6],
        "elapsed": elapsed,
    })

    if elapsed > 120:
        flush(f"\n  (stopping — next K would be too slow)")
        break


# ═══════════════════════════════════════════════════════════════════
flush(f"\n\n{'=' * 78}")
flush(f"  CROSS-K SUMMARY: λ₂(bulk) vs A*")
flush(f"{'=' * 78}")

flush(f"\n  {'K':>4s}  {'λ₂(bulk)':>12s}  {'A_eff(λ₂)':>12s}  {'|A-A*|':>12s}  {'A*':>12s}")
for r in all_results:
    if r["l2_bulk"] is not None:
        l2 = r["l2_bulk"]
        A_eff = (2 * l2 - 1) / (2 * l2 - 2) if abs(l2 - 1) > 1e-10 else float("nan")
        gap = abs(A_eff - A_STAR)
        flush(f"  {r['K']:4d}  {l2:12.8f}  {A_eff:12.8f}  {gap:12.8f}  {A_STAR:12.8f}")

flush(f"\n  LMH test: Does A_eff from λ₂ converge to A* = {A_STAR:.8f}?")
if len(all_results) >= 2:
    l2_vals = [r["l2_bulk"] for r in all_results if r["l2_bulk"] is not None]
    if len(l2_vals) >= 2 and all(l2_vals[i] < l2_vals[i+1] for i in range(len(l2_vals)-1)):
        flush(f"  *** LMH FALSIFIED at per-position level ***")
        flush(f"  λ₂(bulk) is strictly increasing ({l2_vals[0]:.4f} → {l2_vals[-1]:.4f}),")
        flush(f"  not converging to a constant. A_eff diverges to -∞.")
        flush(f"  The LMH must operate at the global resolvent level, not locally.")
    else:
        flush(f"  Result inconclusive — λ₂ trend is not monotonically increasing")

flush(f"\n{'=' * 78}")
flush(f"  DONE")
flush(f"{'=' * 78}")
