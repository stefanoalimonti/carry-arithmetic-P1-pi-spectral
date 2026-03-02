"""
P1-ST: Stopping-time decomposition of R(K) via Cascade Rigidity
=================================================================

Combines Theorem 6 (Cascade Rigidity) with the E187 universal val formula:
  - tau = stopping time = position of first nonzero carry from MSB
  - val = carries[tau-1] - 1  (by CRT, = c_{n+1} - 1 where tau = D-n)
  - carry[tau] = 1 ALWAYS (E187 theorem)

Exact decomposition:
  sigma_{ab} = SUM_n  P(tau=n | sector ab) * E[val | tau=n, sector ab] * N_{ab}

This script:
  1. Enumerates all D-odd sector (0,0) and (1,0) pairs for K=7..13
  2. Records stopping time tau = M for each pair
  3. Computes P(tau | sector) and E[val | tau, sector]
  4. Decomposes R(K) through the stopping-time series
  5. Analyzes sector asymmetry at each depth
  6. Tests for geometric stopping-time distribution
  7. PSLQ on converged constants
"""

import sys
import time
from collections import defaultdict
from math import pi, log, log2

def flush(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def enumerate_pairs(K):
    """Enumerate all D-odd pairs in sectors (0,0) and (1,0).
    Returns per-sector stopping-time data."""
    D = 2 * K - 1
    lo = 1 << (K - 1)
    hi = 1 << K

    # tau_data[sector][(tau, val)] = count
    # sector: "00" or "10"
    tau_counts = {"00": defaultdict(int), "10": defaultdict(int)}
    tau_val_sum = {"00": defaultdict(int), "10": defaultdict(int)}
    sector_sigma = {"00": 0, "10": 0}
    sector_n = {"00": 0, "10": 0}

    for p in range(lo, hi):
        a1 = (p >> (K - 2)) & 1
        for q in range(lo, hi):
            prod = p * q
            if prod.bit_length() != D:
                continue
            b1 = (q >> (K - 2)) & 1

            sec_bits = a1 * 2 + b1
            if sec_bits == 3:  # sector (1,1) is D-odd impossible
                continue
            if sec_bits == 1:  # sector (0,1) — skip, symmetric to (1,0)
                continue

            sec = "10" if sec_bits == 2 else "00"

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
            tau = D - M  # depth from MSB: tau=0 means M=D, tau=2 means M=D-2

            tau_counts[sec][tau] += 1
            tau_val_sum[sec][tau] += val
            sector_sigma[sec] += val
            sector_n[sec] += 1

    return {
        "tau_counts": tau_counts,
        "tau_val_sum": tau_val_sum,
        "sigma": sector_sigma,
        "n": sector_n,
        "D": D,
        "K": K,
    }


# ═══════════════════════════════════════════════════════════════════
flush("=" * 78)
flush("  P1-ST: STOPPING-TIME DECOMPOSITION OF R(K)")
flush("  Cascade Rigidity Theorem + E187 Universal Val Formula")
flush("=" * 78)

all_results = []

for K in range(7, 22):
    D = 2 * K - 1
    lo = 1 << (K - 1)
    hi = 1 << K
    n_total = (hi - lo) ** 2

    flush(f"\n{'─' * 78}")
    flush(f"  K={K}  D={D}  total_pairs={n_total:,}")
    flush(f"{'─' * 78}")

    t0 = time.time()
    data = enumerate_pairs(K)
    elapsed = time.time() - t0
    flush(f"  Enumerated in {elapsed:.1f}s")

    sigma_00 = data["sigma"]["00"]
    sigma_10 = data["sigma"]["10"]
    n_00 = data["n"]["00"]
    n_10 = data["n"]["10"]

    R_K = sigma_10 / sigma_00 if sigma_00 != 0 else float("nan")

    flush(f"\n  Sector totals:")
    flush(f"    σ₀₀ = {sigma_00:>12,d}   n₀₀ = {n_00:>10,d}   <val>₀₀ = {sigma_00/n_00 if n_00 else 0:.6f}")
    flush(f"    σ₁₀ = {sigma_10:>12,d}   n₁₀ = {n_10:>10,d}   <val>₁₀ = {sigma_10/n_10 if n_10 else 0:.6f}")
    flush(f"    R(K) = σ₁₀/σ₀₀ = {R_K:.8f}   (target: -π = {-pi:.8f})")

    # ── Stopping-time distribution ────────────────────────────────
    flush(f"\n  Stopping-time distribution P(τ=n | sector):")
    flush(f"  {'τ':>4s}  {'n₀₀(τ)':>10s}  {'P(τ|00)':>10s}  {'n₁₀(τ)':>10s}  {'P(τ|10)':>10s}"
          f"  {'E[v|τ,00]':>10s}  {'E[v|τ,10]':>10s}  {'ratio_n':>10s}")

    max_tau = max(
        max(data["tau_counts"]["00"].keys(), default=0),
        max(data["tau_counts"]["10"].keys(), default=0),
    )

    tau_data_K = []
    recon_sigma_00 = 0
    recon_sigma_10 = 0

    for tau in range(0, max_tau + 1):
        n00_tau = data["tau_counts"]["00"].get(tau, 0)
        n10_tau = data["tau_counts"]["10"].get(tau, 0)
        vs00 = data["tau_val_sum"]["00"].get(tau, 0)
        vs10 = data["tau_val_sum"]["10"].get(tau, 0)

        p00 = n00_tau / n_00 if n_00 > 0 else 0
        p10 = n10_tau / n_10 if n_10 > 0 else 0
        ev00 = vs00 / n00_tau if n00_tau > 0 else float("nan")
        ev10 = vs10 / n10_tau if n10_tau > 0 else float("nan")
        ratio_n = n10_tau / n00_tau if n00_tau > 0 else float("nan")

        recon_sigma_00 += vs00
        recon_sigma_10 += vs10

        if n00_tau > 0 or n10_tau > 0:
            ev00_s = f"{ev00:10.6f}" if n00_tau > 0 else f"{'—':>10s}"
            ev10_s = f"{ev10:10.6f}" if n10_tau > 0 else f"{'—':>10s}"
            rn_s = f"{ratio_n:10.6f}" if n00_tau > 0 and n10_tau > 0 else f"{'—':>10s}"

            flush(f"  {tau:4d}  {n00_tau:10d}  {p00:10.6f}  {n10_tau:10d}  {p10:10.6f}"
                  f"  {ev00_s}  {ev10_s}  {rn_s}")

        tau_data_K.append({
            "tau": tau, "n00": n00_tau, "n10": n10_tau,
            "p00": p00, "p10": p10, "ev00": ev00, "ev10": ev10,
            "vs00": vs00, "vs10": vs10,
        })

    # Verify reconstruction
    flush(f"\n  Reconstruction check:")
    flush(f"    Σ val_sum(00) = {recon_sigma_00:,d}  (should be {sigma_00:,d})  {'OK' if recon_sigma_00 == sigma_00 else 'MISMATCH!'}")
    flush(f"    Σ val_sum(10) = {recon_sigma_10:,d}  (should be {sigma_10:,d})  {'OK' if recon_sigma_10 == sigma_10 else 'MISMATCH!'}")

    # ── Geometric fit for P(tau|00) ───────────────────────────────
    flush(f"\n  Geometric fit: P(τ=n|00) ~ p₀ * r^n ?")
    p00_vals = [(t["tau"], t["p00"]) for t in tau_data_K if t["n00"] > 10]
    if len(p00_vals) >= 3:
        ratios_00 = []
        for i in range(1, len(p00_vals)):
            if p00_vals[i-1][1] > 1e-10:
                r = p00_vals[i][1] / p00_vals[i-1][1]
                ratios_00.append((p00_vals[i][0], r))
        if ratios_00:
            flush(f"    Successive ratios P(τ=n)/P(τ=n-1) for sector 00:")
            for tau_v, r in ratios_00:
                flush(f"      τ={tau_v:2d}: ratio = {r:.6f}  (1/2 = {0.5:.6f})")

    # ── E[val|tau] convergence across K ───────────────────────────
    flush(f"\n  E[val|τ, sector] by depth:")
    for tau_d in tau_data_K:
        tau = tau_d["tau"]
        if tau_d["n00"] > 10:
            ev00 = tau_d["ev00"]
            crt_check = ""
            if tau >= 2 and tau_d["n00"] > 0:
                crt_check = f"  (CRT: val = c_{tau+1} - 1)"
            flush(f"    τ={tau:2d}  E[val|τ,00]={ev00:+.6f}{crt_check}")

    # ── Sector asymmetry ──────────────────────────────────────────
    flush(f"\n  Sector asymmetry at each depth:")
    flush(f"  {'τ':>4s}  {'P(τ|10)/P(τ|00)':>16s}  {'E[v|10]-E[v|00]':>16s}  {'σ₁₀(τ)/σ₀₀(τ)':>16s}")
    for td in tau_data_K:
        tau = td["tau"]
        if td["n00"] > 10:
            p_ratio = td["p10"] / td["p00"] if td["p00"] > 1e-12 else float("nan")
            ev_diff = (td["ev10"] - td["ev00"]) if td["n10"] > 0 else float("nan")
            s_ratio = td["vs10"] / td["vs00"] if td["vs00"] != 0 else float("nan")

            pr_s = f"{p_ratio:16.6f}" if p_ratio == p_ratio else f"{'—':>16s}"
            evd_s = f"{ev_diff:+16.6f}" if ev_diff == ev_diff else f"{'—':>16s}"
            sr_s = f"{s_ratio:+16.6f}" if s_ratio == s_ratio and abs(s_ratio) < 1e6 else f"{'—':>16s}"
            flush(f"  {tau:4d}  {pr_s}  {evd_s}  {sr_s}")

    all_results.append({
        "K": K, "D": D, "R_K": R_K,
        "sigma_00": sigma_00, "sigma_10": sigma_10,
        "n_00": n_00, "n_10": n_10,
        "tau_data": tau_data_K,
        "elapsed": elapsed,
    })

    if elapsed > 180:
        flush(f"\n  (stopping — next K would be too slow)")
        break


# ═══════════════════════════════════════════════════════════════════
# CROSS-K CONVERGENCE ANALYSIS
# ═══════════════════════════════════════════════════════════════════
flush(f"\n\n{'=' * 78}")
flush(f"  CROSS-K CONVERGENCE ANALYSIS")
flush(f"{'=' * 78}")

# 1. P(tau=n|sector) convergence
flush(f"\n  1. P(τ=n|00) convergence across K:")
flush(f"  {'τ':>4s}", end="")
for r in all_results:
    flush(f"  {'K='+str(r['K']):>10s}", end="")
flush()

for tau in range(0, 10):
    flush(f"  {tau:4d}", end="")
    for r in all_results:
        td = [t for t in r["tau_data"] if t["tau"] == tau]
        if td and td[0]["n00"] > 0:
            flush(f"  {td[0]['p00']:10.6f}", end="")
        else:
            flush(f"  {'—':>10s}", end="")
    flush()

flush(f"\n  2. P(τ=n|10) convergence across K:")
flush(f"  {'τ':>4s}", end="")
for r in all_results:
    flush(f"  {'K='+str(r['K']):>10s}", end="")
flush()

for tau in range(0, 10):
    flush(f"  {tau:4d}", end="")
    for r in all_results:
        td = [t for t in r["tau_data"] if t["tau"] == tau]
        if td and td[0]["n10"] > 0:
            flush(f"  {td[0]['p10']:10.6f}", end="")
        else:
            flush(f"  {'—':>10s}", end="")
    flush()

# 2. E[val|tau,sector] convergence
flush(f"\n  3. E[val|τ,00] convergence across K:")
flush(f"  {'τ':>4s}", end="")
for r in all_results:
    flush(f"  {'K='+str(r['K']):>10s}", end="")
flush()

for tau in range(0, 10):
    flush(f"  {tau:4d}", end="")
    for r in all_results:
        td = [t for t in r["tau_data"] if t["tau"] == tau]
        if td and td[0]["n00"] > 10:
            flush(f"  {td[0]['ev00']:+10.6f}", end="")
        else:
            flush(f"  {'—':>10s}", end="")
    flush()

flush(f"\n  4. E[val|τ,10] convergence across K:")
flush(f"  {'τ':>4s}", end="")
for r in all_results:
    flush(f"  {'K='+str(r['K']):>10s}", end="")
flush()

for tau in range(0, 10):
    flush(f"  {tau:4d}", end="")
    for r in all_results:
        td = [t for t in r["tau_data"] if t["tau"] == tau]
        if td and td[0]["n10"] > 10:
            flush(f"  {td[0]['ev10']:+10.6f}", end="")
        else:
            flush(f"  {'—':>10s}", end="")
    flush()

# 3. Sector asymmetry ratio P(tau|10)/P(tau|00)
flush(f"\n  5. Sector count ratio n₁₀(τ)/n₀₀(τ) convergence:")
flush(f"  {'τ':>4s}", end="")
for r in all_results:
    flush(f"  {'K='+str(r['K']):>10s}", end="")
flush()

for tau in range(0, 10):
    flush(f"  {tau:4d}", end="")
    for r in all_results:
        td = [t for t in r["tau_data"] if t["tau"] == tau]
        if td and td[0]["n00"] > 10 and td[0]["n10"] > 0:
            ratio = td[0]["n10"] / td[0]["n00"]
            flush(f"  {ratio:10.6f}", end="")
        else:
            flush(f"  {'—':>10s}", end="")
    flush()

# 4. Per-depth R contribution
flush(f"\n  6. Per-depth contribution to R(K) = Σ σ₁₀(τ)/σ₀₀:")
flush(f"  {'τ':>4s}", end="")
for r in all_results:
    flush(f"  {'K='+str(r['K']):>10s}", end="")
flush()

for tau in range(0, 10):
    flush(f"  {tau:4d}", end="")
    for r in all_results:
        td = [t for t in r["tau_data"] if t["tau"] == tau]
        if td and r["sigma_00"] != 0:
            contrib = td[0]["vs10"] / r["sigma_00"]
            flush(f"  {contrib:+10.6f}", end="")
        else:
            flush(f"  {'—':>10s}", end="")
    flush()

# 5. Cumulative R reconstruction
flush(f"\n  7. Cumulative R(K) by stopping time:")
flush(f"  {'τ≤':>4s}", end="")
for r in all_results:
    flush(f"  {'K='+str(r['K']):>10s}", end="")
flush()

for tau_max in range(0, 10):
    flush(f"  {tau_max:4d}", end="")
    for r in all_results:
        cum_10 = sum(t["vs10"] for t in r["tau_data"] if t["tau"] <= tau_max)
        cum_00 = sum(t["vs00"] for t in r["tau_data"] if t["tau"] <= tau_max)
        if cum_00 != 0:
            flush(f"  {cum_10/cum_00:+10.4f}", end="")
        else:
            flush(f"  {'—':>10s}", end="")
    flush()

# 6. R(K) summary
flush(f"\n  8. R(K) convergence summary:")
flush(f"  {'K':>4s}  {'R(K)':>12s}  {'|R+π|':>12s}  {'digits':>8s}")
for r in all_results:
    gap = abs(r["R_K"] + pi)
    digits = -log(gap) / log(10) if gap > 0 else float("inf")
    flush(f"  {r['K']:4d}  {r['R_K']:+12.6f}  {gap:12.2e}  {digits:8.2f}")


# Markdown-formatted summary table (for direct paper inclusion)
flush(f"\n  9. Markdown table (§8.2a format, last K):")
r = all_results[-1]
flush(f"| $\\tau$ | $P(\\tau\\mid 00)$ | $P(\\tau\\mid 10)$ | $E[v\\mid\\tau,00]$ | $E[v\\mid\\tau,10]$ | $\\Delta(\\tau)$ |")
flush(f"|--------|----------------:|----------------:|-----------------:|-----------------:|-------------:|")
for tau in range(0, 10):
    td = [t for t in r["tau_data"] if t["tau"] == tau]
    if not td:
        continue
    t = td[0]
    p00 = f"{t['p00']:.3f}" if t['n00'] > 0 else "—"
    p10 = f"{t['p10']:.3f}" if t['n10'] > 0 else "0.000"
    ev00 = f"{t['ev00']:+.3f}" if t['n00'] > 10 else "—"
    ev10 = f"{t['ev10']:+.3f}" if t['n10'] > 10 else "—"
    if t['n00'] > 10 and t['n10'] > 10:
        delta = f"{t['ev10'] - t['ev00']:+.3f}"
    else:
        delta = "—"
    flush(f"| {tau} | {p00} | {p10} | {ev00} | {ev10} | {delta} |")

flush(f"\n{'=' * 78}")
flush(f"  DONE")
flush(f"{'=' * 78}")
