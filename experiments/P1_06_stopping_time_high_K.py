"""
P1-ST-HK: Stopping-time decomposition from E45 high-K data (K=19,20,21)
=========================================================================

Reads the E45 per-position cascade profiles (stop_ab[M], cas_ab[M]) and
converts them to the stopping-time notation of §8.2a:

    tau = D - M    (depth from MSB)
    P(tau | sector) = stop_ab[M] / n_ab
    E[val | tau, sector] = cas_ab[M] / stop_ab[M]

This extends the P1_03 analysis (K=7..14) to K=19,20,21 using the exact
enumeration data from E45_high_K_profiles.c.

Usage:
    python P1_06_stopping_time_high_K.py
"""

import re
import math
import os

PI = math.pi


def parse_e45_file(filename):
    """Parse an E45 output file into a list of K-results."""
    results = []
    current = None

    with open(filename) as f:
        for line in f:
            line = line.rstrip()

            m = re.match(r'>>> COMPLETED K=(\d+)', line)
            if m:
                current = {'K': int(m.group(1)),
                           'cas_00': {}, 'cas_10': {},
                           'stop_00': {}, 'stop_10': {}}
                results.append(current)
                continue

            if current is None:
                continue

            m = re.match(r'\s+D\s+=\s+(\d+)', line)
            if m:
                current['D'] = int(m.group(1))

            m = re.match(r'\s+S00\s+=\s+(-?\d+)\s+\(n00\s*=\s*(\d+)\)', line)
            if m:
                current['S00'] = int(m.group(1))
                current['n00'] = int(m.group(2))

            m = re.match(r'\s+S10\s+=\s+(-?\d+)\s+\(n10\s*=\s*(\d+)\)', line)
            if m:
                current['S10'] = int(m.group(1))
                current['n10'] = int(m.group(2))

            m = re.match(r'\s+R\s+=\s+([+-]?\d+\.\d+)', line)
            if m:
                current['R'] = float(m.group(1))

            m = re.match(r'\s+(\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s*$',
                         line)
            if m:
                d = int(m.group(1))
                current['cas_00'][d] = int(m.group(2))
                current['cas_10'][d] = int(m.group(3))
                current['stop_00'][d] = int(m.group(4))
                current['stop_10'][d] = int(m.group(5))

    return results


def analyze_stopping_time(res):
    """Convert E45 position-based data to stopping-time decomposition."""
    K = res['K']
    D = res['D']
    n00 = res['n00']
    n10 = res['n10']

    tau_data = []
    for tau in range(2, D - 2):
        M = D - tau
        s00 = res['stop_00'].get(M, 0)
        s10 = res['stop_10'].get(M, 0)
        c00 = res['cas_00'].get(M, 0)
        c10 = res['cas_10'].get(M, 0)

        p00 = s00 / n00 if n00 > 0 else 0
        p10 = s10 / n10 if n10 > 0 else 0
        ev00 = c00 / s00 if s00 > 0 else None
        ev10 = c10 / s10 if s10 > 0 else None

        if s00 > 0 or s10 > 0:
            tau_data.append({
                'tau': tau, 'M': M,
                's00': s00, 's10': s10,
                'c00': c00, 'c10': c10,
                'p00': p00, 'p10': p10,
                'ev00': ev00, 'ev10': ev10,
            })

    return tau_data


def main():
    print("=" * 78)
    print("  P1-ST-HK: STOPPING-TIME DECOMPOSITION FROM E45 HIGH-K DATA")
    print("=" * 78)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    e45_dir = os.path.join(os.path.dirname(os.path.dirname(script_dir)),
                           "carry-arithmetic-E-trace-anomaly", "experiments")

    all_results = []
    for fname in ["E45_K19_K20.txt", "E45_K21.txt"]:
        fpath = os.path.join(e45_dir, fname)
        if os.path.exists(fpath):
            all_results.extend(parse_e45_file(fpath))
            print(f"  Loaded {fname}")
        else:
            print(f"  WARNING: {fname} not found at {fpath}")

    if not all_results:
        print("  No E45 data found. Exiting.")
        return

    for res in all_results:
        K = res['K']
        D = res['D']
        R = res.get('R', res['S10'] / res['S00'] if res['S00'] != 0 else float('nan'))
        n00 = res['n00']
        n10 = res['n10']

        print(f"\n{'─' * 78}")
        print(f"  K = {K},  D = {D},  R = {R:+.12f},  |R+π| = {abs(R + PI):.6e}")
        print(f"  n00 = {n00:,},  n10 = {n10:,}")

        tau_data = analyze_stopping_time(res)

        # Summary table for small tau (boundary layer)
        print(f"\n  Stopping-time decomposition (τ = D − M, boundary layer):")
        print(f"  {'τ':>4s}  {'P(τ|00)':>10s}  {'P(τ|10)':>10s}  "
              f"{'E[v|τ,00]':>10s}  {'E[v|τ,10]':>10s}  {'Δ(τ)':>10s}  "
              f"{'n00':>10s}  {'n10':>10s}")

        for td in tau_data:
            tau = td['tau']
            if tau > 12:
                continue

            p00_s = f"{td['p00']:.6f}" if td['s00'] > 0 else "—"
            p10_s = f"{td['p10']:.6f}" if td['s10'] > 0 else "0.000000"
            ev00_s = f"{td['ev00']:+.6f}" if td['ev00'] is not None and td['s00'] > 100 else "—"
            ev10_s = f"{td['ev10']:+.6f}" if td['ev10'] is not None and td['s10'] > 100 else "—"

            if td['ev00'] is not None and td['ev10'] is not None and td['s00'] > 100 and td['s10'] > 100:
                delta = td['ev10'] - td['ev00']
                delta_s = f"{delta:+.6f}"
            else:
                delta_s = "—"

            print(f"  {tau:4d}  {p00_s:>10s}  {p10_s:>10s}  "
                  f"{ev00_s:>10s}  {ev10_s:>10s}  {delta_s:>10s}  "
                  f"{td['s00']:>10,}  {td['s10']:>10,}")

        # Cumulative R reconstruction
        print(f"\n  Cumulative R through stopping depth:")
        cum_c10, cum_c00 = 0, 0
        for td in tau_data:
            cum_c10 += td['c10']
            cum_c00 += td['c00']
            if td['tau'] <= 12 and cum_c00 != 0:
                R_cum = cum_c10 / cum_c00
                print(f"    τ ≤ {td['tau']:2d}:  R_cum = {R_cum:+.6f}  "
                      f"(|R_cum + π| = {abs(R_cum + PI):.4e})")

        # Full R check
        total_c00 = sum(td['c00'] for td in tau_data)
        total_c10 = sum(td['c10'] for td in tau_data)
        if total_c00 != 0:
            R_full = total_c10 / total_c00
            print(f"    FULL:     R_full = {R_full:+.12f}  "
                  f"(|R + π| = {abs(R_full + PI):.6e})")
            print(f"    E45 R:    R      = {R:+.12f}")

    # Cross-K comparison of converged constants
    print(f"\n{'=' * 78}")
    print(f"  CROSS-K CONVERGENCE OF UNIVERSAL CONSTANTS")
    print(f"{'=' * 78}")
    print(f"\n  {'τ':>4s}", end="")
    for res in all_results:
        print(f"  {'K='+str(res['K']):>14s}", end="")
    print()

    print(f"\n  --- P(τ|00) ---")
    for tau in range(2, 10):
        print(f"  {tau:4d}", end="")
        for res in all_results:
            td_list = analyze_stopping_time(res)
            td = [t for t in td_list if t['tau'] == tau]
            if td and td[0]['s00'] > 0:
                print(f"  {td[0]['p00']:14.8f}", end="")
            else:
                print(f"  {'—':>14s}", end="")
        print()

    print(f"\n  --- P(τ|10) ---")
    for tau in range(2, 10):
        print(f"  {tau:4d}", end="")
        for res in all_results:
            td_list = analyze_stopping_time(res)
            td = [t for t in td_list if t['tau'] == tau]
            if td and td[0]['s10'] > 0:
                print(f"  {td[0]['p10']:14.8f}", end="")
            else:
                print(f"  {'—':>14s}", end="")
        print()

    print(f"\n  --- E[val|τ,00] ---")
    for tau in range(2, 10):
        print(f"  {tau:4d}", end="")
        for res in all_results:
            td_list = analyze_stopping_time(res)
            td = [t for t in td_list if t['tau'] == tau]
            if td and td[0]['ev00'] is not None and td[0]['s00'] > 100:
                print(f"  {td[0]['ev00']:+14.8f}", end="")
            else:
                print(f"  {'—':>14s}", end="")
        print()

    print(f"\n  --- E[val|τ,10] ---")
    for tau in range(2, 10):
        print(f"  {tau:4d}", end="")
        for res in all_results:
            td_list = analyze_stopping_time(res)
            td = [t for t in td_list if t['tau'] == tau]
            if td and td[0]['ev10'] is not None and td[0]['s10'] > 100:
                print(f"  {td[0]['ev10']:+14.8f}", end="")
            else:
                print(f"  {'—':>14s}", end="")
        print()

    print(f"\n  --- Δ(τ) = E[v|τ,10] - E[v|τ,00] ---")
    for tau in range(2, 10):
        print(f"  {tau:4d}", end="")
        for res in all_results:
            td_list = analyze_stopping_time(res)
            td = [t for t in td_list if t['tau'] == tau]
            if (td and td[0]['ev00'] is not None and td[0]['ev10'] is not None
                    and td[0]['s00'] > 100 and td[0]['s10'] > 100):
                delta = td[0]['ev10'] - td[0]['ev00']
                print(f"  {delta:+14.8f}", end="")
            else:
                print(f"  {'—':>14s}", end="")
        print()

    # Markdown table for paper
    print(f"\n{'=' * 78}")
    print(f"  MARKDOWN TABLE (§8.2a format, K=21 data)")
    print(f"{'=' * 78}")
    res21 = [r for r in all_results if r['K'] == 21]
    if res21:
        res = res21[0]
        td_list = analyze_stopping_time(res)
        print(f"| $\\tau$ | $P(\\tau\\mid 00)$ | $P(\\tau\\mid 10)$ | "
              f"$E[v\\mid\\tau,00]$ | $E[v\\mid\\tau,10]$ | $\\Delta(\\tau)$ |")
        print(f"|--------|----------------:|----------------:|"
              f"-----------------:|-----------------:|-------------:|")
        for td in td_list:
            if td['tau'] > 10:
                break
            p00 = f"{td['p00']:.3f}" if td['s00'] > 0 else "—"
            p10 = f"{td['p10']:.3f}" if td['s10'] > 0 else "0.000"
            ev00 = f"{td['ev00']:+.3f}" if td['ev00'] is not None and td['s00'] > 100 else "—"
            ev10 = f"{td['ev10']:+.3f}" if td['ev10'] is not None and td['s10'] > 100 else "—"
            if (td['ev00'] is not None and td['ev10'] is not None
                    and td['s00'] > 100 and td['s10'] > 100):
                delta = f"{td['ev10'] - td['ev00']:+.3f}"
            else:
                delta = "—"
            print(f"| {td['tau']} | {p00} | {p10} | {ev00} | {ev10} | {delta} |")

    print(f"\n{'=' * 78}")
    print(f"  DONE")
    print(f"{'=' * 78}")


if __name__ == '__main__':
    main()
