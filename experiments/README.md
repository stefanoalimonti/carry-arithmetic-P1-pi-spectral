# Experiments — Paper P1

| Script | Description | Referenced in |
|--------|-------------|---------------|
| `P1_01_bernoulli_hierarchy.py` | Bernoulli hierarchy decomposition by carries at D−2, D−3, D−4 (K=7..13). Confirms B₃=0 and discovers cascade value rigidity. | §8.1 (Levels 3–4) |
| `P1_02_cascade_rigidity.py` | Computational verification of Theorem 6 (Cascade Rigidity): checks val determinism in all (c₂,c₃,c₄) sub-cells, K=7..13. | §8.1 (Theorem 6) |
| `P1_03_stopping_time_decomposition.py` | Stopping-time decomposition of R(K): records τ, P(τ\|sector), E[val\|τ,sector], sector asymmetry Δ(τ). Outputs markdown table. K=7..14. | §8.2a |
| `P1_05_doob_eigenvalues.py` | Doob h-transform eigenvalues: per-position carry transition matrices, bulk eigenvalues. **Falsifies LMH at per-position level** (λ₂ grows 0.69→0.86). K=5..12. | §8.2b |
| `P1_06_stopping_time_high_K.py` | Stopping-time decomposition from E45 high-K data (K=19,20,21). Extends P1\_03 to 5-digit precision; PSLQ on extrapolated constants; markdown table output. | §8.2a |
| `P1_07_richardson_extrapolation.py` | Exponential Richardson extrapolation of R(K): Neville–Aitken table with h=(1/2)^K, Aitken Δ², leave-one-out cross-validation. Reproduces the 4.0–4.4-digit convergence claim (depending on input precision). | §9.2 (Proposition 3) |

Additional computational verification for P1 is in the companion repository [`carry-arithmetic-E-trace-anomaly`](https://github.com/stefanoalimonti/carry-arithmetic-E-trace-anomaly), experiments E09–E30.

## Requirements

Python >= 3.8, NumPy.
