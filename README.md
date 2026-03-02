# carry-arithmetic-P1-pi-spectral

**Pi from Pure Arithmetic: A Spectral Phase Transition in the Binary Carry Bridge**

*Author: Stefano Alimonti* · [ORCID 0009-0009-1183-1698](https://orcid.org/0009-0009-1183-1698)

## Main Result

Numerically, the sector-ratio trace of binary multiplication converges toward $-\pi$ (4.0-digit exponential Richardson extrapolation through K = 21). The paper develops five equivalent formulations of this conjecture ($R = -\pi$, Conjecture 1): series identity over depths, angular integral, carry weight computation, resolvent trace, and separable weight matrix with carry constraints. Theorems 1–7 are proved unconditionally; Conjecture 1 ($R = -\pi$) remains open.

**Theorem 6 (Cascade Rigidity):** In sector (0,0), the cascade value is determined algebraically by the near-boundary carry profile: if $c_n \geq 1$ while $c_2 = \cdots = c_{n-1} = 0$, then $M = D-n$ and $\text{val} = c_{n+1} - 1$. This proves $B_3 = 0$ in the Bernoulli hierarchy and extends recursively to arbitrary depth.

**Stopping-time decomposition (§8.2a):** Combining Theorem 6 with the unit leading carry, R(K) decomposes exactly as a sum over stopping depths τ. Both $P(\tau \mid \text{sector})$ and $E[\text{val} \mid \tau, \text{sector}]$ converge to K-independent universal constants. The sector asymmetry $\Delta(\tau)$ forms a stable sequence (~0.387, 0.351, 0.411, 0.315, 0.335). ~95% of R(K) concentrates in τ ≤ 7.

## Status

- **Proved:** structural framework, Theorems 1–7, exact stopping-time resolvent identities.
- **Conditional:** closed-form `R = S(A)` under LMH.
- **Open target:** scalar closure `C=-4` in `R(∞)=C·L(1,χ₄)` and analytic proof of dominant `1/2` rate preservation for the conditioned chain.

Main paper (~20pp). The carry-Dirichlet channel extending the stopping-time framework to a function of s is developed in [L].

## Repository Structure

```
paper/pi_from_arithmetic.md                      The paper
experiments/
  P1_01_bernoulli_hierarchy.py                   Bernoulli hierarchy decomposition (§8.1)
  P1_02_cascade_rigidity.py                      Cascade Rigidity Theorem verification (§8.1)
  P1_03_stopping_time_decomposition.py           Stopping-time decomposition of R(K) (§8.2a)
  P1_05_doob_eigenvalues.py                      Doob h-transform eigenvalue analysis (§8.2b)
  P1_06_stopping_time_high_K.py                  High-K stopping-time from E45 data (§8.2a)
```

Experiments P1_01–P1_03 verify Theorem 6, the Bernoulli hierarchy, and the stopping-time decomposition (§8.1–§8.2a). P1_05 tests — and falsifies — the LMH at the per-position eigenvalue level (§8.2b). Additional computational verification for Conjectures 1–2 and Theorems 1–3, 5 is in the companion repository `carry-arithmetic-E-trace-anomaly`, experiments E09–E30.

## Reproduction

```bash
pip install numpy
python experiments/P1_01_bernoulli_hierarchy.py   # B₃=0 confirmation, cascade rigidity
python experiments/P1_02_cascade_rigidity.py       # Theorem 6 verification (K=7..13)
python experiments/P1_03_stopping_time_decomposition.py  # Stopping-time decomposition (K=7..14)
python experiments/P1_05_doob_eigenvalues.py       # Doob eigenvalues — LMH falsified (K=5..12)
```

## Dependencies

- Python >= 3.8, NumPy

## Companion Papers

| Label | Title | Repository |
|-------|-------|------------|
| [E] | The Trace Anomaly of Binary Multiplication (**experiments**) | [`carry-arithmetic-E-trace-anomaly`](https://github.com/stefanoalimonti/carry-arithmetic-E-trace-anomaly) |
| [P2] | The Sector Ratio in Binary Multiplication | [`carry-arithmetic-P2-sector-ratio`](https://github.com/stefanoalimonti/carry-arithmetic-P2-sector-ratio) |
| [G] | The Angular Uniqueness of Base 2 | [`carry-arithmetic-G-angular-uniqueness`](https://github.com/stefanoalimonti/carry-arithmetic-G-angular-uniqueness) |
| [L] | The Carry–Dirichlet Bridge | [`carry-arithmetic-L-dirichlet-bridge`](https://github.com/stefanoalimonti/carry-arithmetic-L-dirichlet-bridge) |

### Citation

```bibtex
@article{alimonti2026pi_arithmetic,
  author  = {Alimonti, Stefano},
  title   = {Pi from Pure Arithmetic: A Spectral Phase Transition in the Binary Carry Bridge},
  year    = {2026},
  note    = {Preprint},
  url     = {https://github.com/stefanoalimonti/carry-arithmetic-P1-pi-spectral}
}
```

## License

Paper: CC BY 4.0. Code: MIT License.
