# π from Pure Arithmetic: A Spectral Phase Transition in the Binary Carry Bridge

**Author:** Stefano Alimonti
**Affiliation:** Independent Researcher
**Date:** March 2026

---

## Abstract

We study the propagation of arithmetic carries in binary multiplication and discover a spectral phase transition that converts rational arithmetic into a transcendental constant. Given two uniformly random K-bit integers X, Y with product P = XY written in D = 2K − 1 bits, define the *sector ratio* R(K) = σ₁₀(K)/σ₀₀(K), where σ_ac is the carry-weighted partition function restricted to pairs whose second-highest bits of X, Y are (a, c).

We prove:

1. The Markov (independent-convolution) model gives R_Markov = +2/3, a rational constant independent of K (Theorem 2).
2. The alternating spectral sum over the Dirichlet bridge has the exact closed form S(β) = 2β/(1 + 2β), a Möbius transformation (Theorem 1).
3. Setting S = −π yields the critical parameter A* = (3π + 2)/(2(π + 1)), placing the system at distance 1/(2(π + 1)) from the pole (Corollary 2).
4. The spectral zeta function of the Markov carry bridge reproduces ζ(2s) at all even arguments of the Riemann zeta function, via Weyl asymptotics (Theorem 3).
5. No D-odd pair has sector (1,1) (Theorem 4), and the count ratio N₁₀/N₀₀ → (2ln(4/3) − 1/2)/(2ln(9/8)) = 0.31993… involves only ln 2 and ln 3 (Theorem 5). The entire factor π in R enters exclusively through the carry weight ratio ⟨val⟩₁₀/⟨val⟩₀₀.

We conjecture:

6. The exact cascade sector ratio satisfies R(∞) = lim_{K→∞} R(K) = −π = −4 L(1, χ₄) (Conjecture 1). Exact enumeration up to K = 21 (1.1 × 10¹² digit pairs) confirms this to 4.0 significant digits by exponential Richardson extrapolation (§9.2; 4.4 digits using exact values only). The weight sin(nπ/2) in the spectral sum equals the Dirichlet character χ₄(n), which emerges as the effective spectral weight under resolvent universality [E].

The system undergoes a spectral phase transition: the Markov model (subcritical) produces the rational +2/3; digit correlations push the carry bridge past a critical point into a supercritical regime where the transcendental −π emerges. Combined with ζ(2s) from the spectral zeta function and L(1, χ₄) from the sector ratio, the carry bridge produces both factors of the Dedekind zeta function ζ_{ℚ(i)}(s) = ζ(s) · L(s, χ₄) of the Gaussian integers.

---

## 1. Introduction

The appearance of π in number theory is usually traced to analytic methods: the residue calculus in the proof of the prime number theorem, the Fourier analysis underlying Dirichlet L-functions, or the Gaussian integral in the central limit theorem. In each case, the geometric content of π (circumference, area, rotation) enters through an analytical tool applied to a discrete problem.

In this paper we exhibit a mechanism that produces π from *purely arithmetic operations* — binary multiplication and carry propagation — without invoking any geometric or analytic machinery except at the final step of evaluating a closed-form identity. The mechanism reveals a spectral phase transition: a rational quantity (+2/3) is converted into a transcendental one (−π) when digit correlations push the system past a critical point.

We also show that the Markov transfer operator for the carry bridge, viewed as a spectral zeta function, reproduces ζ(2s) at all even arguments — a manifestation of Weyl's law in the arithmetic setting. Combined with the conjectural L(1, χ₄) = π/4 from the sector ratio, this connects the carry bridge to the Dedekind zeta function of the Gaussian integers.

**Scope and status.** Five theorems are proved unconditionally (Theorems 1–5, Theorem 6, Theorem 7). The central claim R(∞) = −π (Conjecture 1) is supported by 4.0-digit numerical evidence (exponential Richardson, §9.2) but remains open; under the Linear Mix Hypothesis of [E], it follows from the proved closed-form identity S(A) = 2(1−A)/(3−2A). This paper is therefore a contribution to experimental mathematics: it identifies the phenomenon, proves the structural framework, and isolates the single remaining spectral conjecture. The remaining proof targets are: (i) prove the scalar identity C = −4 in R(∞) = C · L(1, χ₄), and (ii) prove analytically that the conditioned D-odd chain preserves the dominant 1/2 spectral rate. The stopping-time framework extends to a full Dirichlet series in the complex variable s in [L], revealing L²(s, χ₄) structure.

### 1.1 Setting

Let X and Y be independent uniformly random integers in [2^{K−1}, 2^K) (K-bit integers with leading bit 1). Their product P = XY has D = 2K − 1 binary digits. The schoolbook multiplication algorithm computes

    P = Σ_{j=0}^{D-1} p_j · 2^j

where at each position j the product digit and carry are determined by

    p_j = (conv_j + c_j) mod 2,    c_{j+1} = ⌊(conv_j + c_j) / 2⌋

with conv_j = Σ_i x_i · y_{j−i} the convolution at position j and c_0 = 0 the boundary condition.

**Definition 1 (Sector partition).** The *sector* of a pair (X, Y) is (a, c) ∈ {0,1}² where a = x_{K−2} and c = y_{K−2} are the second-highest bits. The *sector partition function* is

    σ_ac(K) = Σ_{(X,Y): x_{K-2}=a, y_{K-2}=c} w(X, Y)

where $w(X,Y) = w_{\mathrm{casc}}(X,Y)$ is the cascade carry weight defined below.  Two valuations appear in the literature: the *schoolbook* weight $w_{\mathrm{sch}} = 2c_{D-2} - 1$ (fixed-position readout) and the *cascade* weight $w_{\mathrm{casc}} = c_{M-1} - 1$ where $M$ is the topmost nonzero carry position (first-passage readout; not to be confused with the companion matrix $M$ of [A], [B]).  Throughout this paper $w = w_{\mathrm{casc}}$ unless otherwise noted; it is this valuation that produces $-\pi$.

**Definition 2 (Sector ratio).** R(K) = σ₁₀(K)/σ₀₀(K). By symmetry, σ₁₀ = σ₀₁.  (Not to be confused with the carry correction factor $R(l,s)$ of [B], which measures the per-prime deviation from the Euler product.)

### 1.2 Main results

**Theorem 1 (Spectral closed form).** For β ∈ (−1/2, +∞), define

    S(β, L) = Σ_{k=0}^{⌊L/2⌋-1} (-1)^k / (1/2 + β sin²((2k+1)π/(2L)))

Then

    lim_{L→∞} S(β, L) = 2β/(1 + 2β)

**Corollary 1 (Möbius form).** Setting β = 1 − A, the sum as a function of the eigenvalue shift A is

    S(A) = 2(1−A)/(3−2A)

a Möbius transformation with a zero at A = 1 and a pole at A = 3/2.

**Corollary 2 (Phase transition).**

- A = 0 (Markov): S = +2/3
- A = 1 (critical): S = 0
- A = A*: S = −π, where A* = (3π + 2)/(2(π + 1)) = 3/2 − 1/(2(π+1)) ≈ 1.37927

**Theorem 2 (Markov baseline).**  The independent-convolution (Markov) model gives R_Markov = +2/3, with finite-size correction

    S(0, L) = 2/3 + 10 ζ(2)/(3 L²) + O(L⁻⁴)

**Theorem 3 (Spectral zeta function).** Let F = ½I − T_bridge be the fluctuation operator of the Markov carry bridge of length L, with eigenvalues ε_n = sin²(nπ/(2(L+1))). Then for Re(s) > 1/2,

    ζ_F(s) / (2L/π)^{2s}  →  ζ(2s)    as L → ∞

with convergence rate $O(L^{1-2\sigma})$ where $\sigma = \Re(s)$ (reducing to $O(1/L^2)$ for $\sigma > 3/2$).

**Theorem 4 (N₁₁ = 0).** No D-odd pair has sector (1,1).

**Theorem 5 (Closed-form count ratio).** N₁₀/N₀₀ → (2ln(4/3) − 1/2)/(2ln(9/8)) = 0.31993…, involving only ln 2 and ln 3. The convergence correction is O(1/2^K).

**Conjecture 1 (Exact sector ratio).** R(∞) = lim_{K→∞} σ₁₀(K)/σ₀₀(K) = −π. (Note: this is implied by, but weaker than, the Linear Mix Hypothesis of [E], which is also called "Conjecture 1" in that paper.)

Supported by exact enumeration to K = 21 with Richardson extrapolation (Section 9). By the decomposition R = R₀ + ΔR [E], the conjecture reduces to proving that the cascade correction ΔR = Tr[(I−T)⁻¹·P_{c=0}] equals exactly −π − R₀ = +0.7896...

### 1.3 Structure of the paper

- Section 2: Carry bridge and Dirichlet eigenfunctions
- Section 3: Markov baseline (R = +2/3)
- Section 4: Proof of the closed form (Theorem 1)
- Section 5: The spectral phase transition
- Section 6: The Dirichlet character χ₄ and L(1, χ₄)
- Section 7: The spectral zeta function and ζ(2s) (Theorem 3)
- Section 8: The Bernoulli hierarchy and the Gaussian integers
- Section 9: Numerical evidence
- Section 10: The gap between S(A) and R(K)
- Section 11: Open problems

---

## 2. The Carry Bridge and Dirichlet Eigenfunctions

### 2.1 The D-odd constraint

A pair (X, Y) with X, Y ∈ [2^{K−1}, 2^K) has product P with exactly D = 2K − 1 bits (the "D-odd" constraint). This forces c_0 = 0 (no carry into position 0) and c_{D−1} = 0 (no overflow beyond D − 1 bits). The carry chain c_0, c_1, …, c_{D−1} with these boundary conditions forms a *carry bridge* — a lattice path conditioned to return to zero, analogous to a Brownian bridge.

### 2.2 Transfer operator and eigenvalues

The Markov transfer operator for the carry chain, conditioned on independent convolutions at each position, was introduced by Diaconis and Fulman [1]. In the base-2 case, the carry at each position takes values in {0, 1}, and the transition matrix has eigenvalues

    λ_n = (1/2) cos(nπ/L),    n = 1, …, L−1

where L = D = 2K − 1 is the bridge length. The corresponding eigenfunctions are the Dirichlet sine functions φ_n(j) = sin(nπj/L), reflecting the zero boundary conditions c_0 = c_L = 0. The leading eigenvalue λ_1 = (1/2)cos(π/L) → 1/2 encodes the Diaconis–Fulman spectral gap.

### 2.3 Spectral representation of the sector response

The sector perturbation at position j* = K − 2 excites the Dirichlet modes of the carry bridge. The response at the boundary is a spectral sum

    R ~ Σ_{n=1}^{L-1} sin(nπj*/L) / (1 − λ_n^eff)

where λ_n^eff incorporates both the Markov eigenvalue and the effective modification due to digit correlations. In the Markov model, the perturbation position satisfies j*/L = (K−2)/(2K−1) → 1/2 as K → ∞, giving sin(nπj*/L) → sin(nπ/2) = χ₄(n), the non-principal Dirichlet character modulo 4 (§6). In the true cascade, the stopping position is not localized at the midpoint (it concentrates near the top of the chain; see [E, §8.5]). Nonetheless, under the LMH, the resolvent-weighted spectral sum converges to −π through a collective mechanism [E, §8.3] that is functionally equivalent to the midpoint evaluation (§6.2).

The precise relationship between R(K) and the spectral sum S(A, L) is subtle: [E, §8.1–§8.3] establish that S(A, L) is not a direct model of R(K), but rather that both quantities converge (numerically) toward −π through different mechanisms (§10). Exact enumeration through K = 21 ($5.5 \times 10^{11}$ pairs; [E, Tables 1–2]) confirms convergence to −π with 2.1 significant digits raw and 4.0 digits by exponential Richardson extrapolation (§9). [E] identifies the mechanism as **resolvent universality** of the effective transfer operator. Formalizing this as a rigorous identity remains the principal open problem.

---

## 3. The Markov Baseline: R = +2/3

### 3.1 Independent convolutions

In the Markov model, convolutions at distinct positions are treated as independent random variables given the carry state. The sector perturbation modifies the convolution distribution at position K − 2 by conditioning on the second-highest bits (a, c) of the input.

### 3.2 Proof of S(0, L) → 2/3

With A = 0 (no correlation-induced shift), the spectral sum reduces to

    S(0, L) = Σ_{n odd} (-1)^{(n-1)/2} / (1 − (1/2)cos(nπ/L))

This is S(β = 1, L) in the notation of Theorem 1, giving S(β = 1) = 2·1/(1+2) = 2/3.

**Proposition 1 (Finite-size correction).**

    S(0, L) − 2/3 = 5π²/(9L²) + O(L⁻⁴) = 10 ζ(2)/(3L²) + O(L⁻⁴)

The correction coefficient 5π²/9 is derived analytically from the Fourier expansion (Section 4) and verified numerically to 6 significant digits.

---

## 4. Proof of the Closed Form

**Proof of Theorem 1.** The proof proceeds in five steps.

**Step 1: Fourier expansion.** The integrand has the standard Fourier series (see e.g. [10])

    1/(1/2 + β sin²t) = (2/γ) [1 + 2 Σ_{m=1}^∞ r^m cos(2mt)]

where γ = √(1 + 2β) and r = (γ − 1)²/(2β), valid for |r| < 1 (equivalently β > −1/2). This follows from the identity for 1/(A − B cos θ) with A = (1 + β), B = β, θ = 2t.

**Step 2: Alternating midpoint sum.** For L ≡ 0 (mod 4), the alternating midpoint sum satisfies

    T_m(L) := Σ_{k=0}^{L/2-1} (-1)^k cos(m(2k+1)π/L) =
        1/cos(mπ/L)   if m is odd,
        0              if m is even.

*Proof.* By geometric summation, Σ_k (-1)^k e^{im(2k+1)π/L} = e^{imπ/L}(1 − (−1)^{L/2+m}) / (1 + e^{2imπ/L}). When L/2 is even and m is odd, the numerator is 2 and the denominator is 2cos(mπ/L). ∎

**Step 3: Assembly.** Substituting the Fourier expansion into S(β, L):

    S(β, L) = (2/γ) [T_0 + 2 Σ_{m=1}^∞ r^m T_m(L)]
            = (4/γ) Σ_{j=0}^∞ r^{2j+1} / cos((2j+1)π/L)

(T_0 = 0 because the alternating sum of a constant vanishes when L/2 is even.)

**Step 4: Continuum limit.** As L → ∞, cos((2j+1)π/L) → 1 for each fixed j. Since |r| < 1, dominated convergence gives

    S(β, ∞) = (4r/γ) · 1/(1 − r²)

**Step 5: Simplification.** Direct computation using r = (γ − 1)²/(2β) and γ² = 1 + 2β:

    4r / (γ(1 − r²)) = 4·(γ−1)²/(2β) / (γ · 4γ(γ−1)/(2β)²)

After cancellation, this simplifies to

    S(β, ∞) = 2β/(1 + 2β)   ∎

**Proposition 2 (General correction coefficient).** The O(1/L²) correction to S(β, L) is

    S(β, L) = 2β/(1 + 2β) + C(β)/L² + O(L⁻⁴)

where

    C(β) = 2π² r(1 + 6r² + r⁴) / (γ(1 − r²)³)

For β = 1: C(1) = 5π²/9 = 10ζ(2)/3, confirming Proposition 1.

*Proof.* Expand 1/cos((2j+1)π/L) = 1 + (2j+1)²π²/(2L²) + O(L⁻⁴), substitute into Step 3, and sum the geometric series. ∎

**Remark (M-parity subtlety).** The closed form S(β,∞) = 2β/(1+2β) applies when M = ⌊L/2⌋ is even. For M odd, the limit differs by an O(1) term at finite L, but both subsequences converge to the same limit. This is because the additional term vanishes as L → ∞.

---

## 5. The Spectral Phase Transition

### 5.1 Interpretation

The closed form S(A) = 2(1−A)/(3−2A) reveals a *spectral phase transition* parameterized by the eigenvalue shift A:

| A | S(A) | Regime |
|---|------|--------|
| 0 | +2/3 | Subcritical (Markov): rational |
| 1/2 | +1/2 | Subcritical |
| 1 | 0 | Critical: exact cancellation |
| 1.2 | −2/3 | Supercritical |
| A* ≈ 1.379 | −π | Supercritical: transcendental |
| 3/2 | ±∞ | Pole (divergence) |

### 5.2 The critical point

At A = 1 (β = 0), all effective eigenvalues equal 1/2, and every term in the alternating sum has the same magnitude. The sum vanishes by exact cancellation: Σ (−1)^k · 2 = 0. This is the boundary between the subcritical (positive) and supercritical (negative) regimes.

### 5.3 Subcritical vs. supercritical

For A < 1, the effective eigenvalues are below 1/2 and the alternating sum is positive (dominated by the k = 0 term). For A > 1, eigenvalues exceed 1/2, denominators shrink, and the negative terms dominate. The sign flip at A = 1 is the *trace anomaly*: digit correlations reverse the sign of the classical (Markov) result.

### 5.4 Distance to the pole

The value A* = 3/2 − 1/(2(π+1)) places the physical system at distance 1/(2(π+1)) ≈ 0.121 from the pole of the Möbius transformation. The formula is self-referential: A* contains π, and S(A*) = −π. This is not circular — it is a well-defined algebraic relation 2(1−A)/(3−2A) = −π with a unique solution A = (3π+2)/(2(π+1)).

---

## 6. The Dirichlet Character and L-function

### 6.1 sin(nπ/2) = χ₄(n)

The weight function in the spectral sum takes the values

    sin(nπ/2) = +1 if n ≡ 1 (mod 4),
               = −1 if n ≡ 3 (mod 4),
               =  0 if n is even,

which is exactly the non-principal Dirichlet character χ₄: (ℤ/4ℤ)* → {±1}.

### 6.2 The effective χ₄ model

In the Markov (independent-convolution) model, the sector perturbation acts at position j* = K − 2, and the bridge has length L = 2K − 1. Since j*/L → 1/2 as K → ∞, the Dirichlet mode evaluation sin(nπj*/L) → sin(nπ/2) = χ₄(n). This "midpoint selection" provides a clean motivation for why χ₄ appears in the spectral sum.

*Caveat ([E]).* In the true (non-Markov) cascade, the stopping position concentrates near the top of the chain (depth j = 2, 3, 4 from position D−2), not at the midpoint. The physical cascade stopping profile, projected onto the Dirichlet sine basis, does not equal χ₄(n) at finite K [E, §8]. Nonetheless, the resolvent

$$(I - \mathcal{K}_{\mathrm{eff}})^{-1}$$

acts as a matched spectral filter that transmutes the true erratic cascade profiles into a macroscopic sum converging (under LMH) to −π [E, §8.3]. The χ₄ model is therefore best understood as an *effective toy reward* — a mathematically tractable proxy that belongs to the same universality class as the true cascade in the thermodynamic limit.

### 6.3 Connection to L(1, χ₄)

The Dirichlet L-function at s = 1 gives the Leibniz formula:

    L(1, χ₄) = Σ_{n=0}^∞ (−1)^n/(2n+1) = 1 − 1/3 + 1/5 − ⋯ = π/4

If R = −π (Conjecture 1), then R = −4 L(1, χ₄). The carry bridge produces π through the same Dirichlet series that appears in the classical evaluation of π via the Leibniz–Gregory formula.

### 6.4 Structural identity

Combined with the sub-sector ratio α → 1/6 = B₂ (the second Bernoulli number) and ζ(2) = π²/6:

    α · R² = (1/6) · π² = π²/6 = ζ(2)

If both α = 1/6 (Conjecture 2, §8) and R = −π (Conjecture 1) can be established from carry arithmetic, this would yield a new proof of the Basel problem — one that begins with integer multiplication. **Caveat:** this depends on two currently unproved conjectures.

---

## 7. The Spectral Zeta Function and ζ(2s)

### 7.1 The fluctuation operator

**Definition 3 (Fluctuation operator).** F = ½I − T_bridge, where T_bridge is the Markov carry bridge transfer matrix. The eigenvalues of F are

    ε_n = 1/2 − λ_n = sin²(nπ/(2(L+1))),    n = 1, …, L−1

**Definition 4 (Spectral zeta function).**

    ζ_F(s) = Σ_{n=1}^{L-1} ε_n^{-s} = Σ_{n=1}^{L-1} sin^{-2s}(nπ/(2(L+1)))

### 7.2 Proof of Theorem 3

For n ≪ L, the small-angle approximation sin(x) ≈ x gives ε_n ≈ (nπ/(2L))², so

    ζ_F(s) ≈ (2L/π)^{2s} Σ_{n=1}^{L-1} n^{-2s}

We make this precise. Write ε_n = (nπ/(2L))² · (sin(x_n)/x_n)² where x_n = nπ/(2L). Then

    ε_n^{-s} = (2L/(nπ))^{2s} · (x_n/sin(x_n))^{2s}

and

    ζ_F(s)/(2L/π)^{2s} = Σ_{n=1}^{L-1} n^{-2s} · (x_n/sin(x_n))^{2s}

**Error analysis.** The error has two sources:

*(i) Tail truncation.* The missing terms Σ_{n≥L} n^{-2s} contribute O(L^{1−2s}) for Re(s) > 1/2.

*(ii) Sine correction.* For each mode, (x/sin(x))^{2s} = 1 + sx²/3 + O(x⁴). Summing the leading correction:

    Σ_{n=1}^{L-1} n^{-2s} · s·(nπ/(2L))²/3 = (sπ²/(12L²)) · Σ_{n=1}^{L-1} n^{2−2s}

For s > 3/2 this sum converges and the correction is O(1/L²). For 1/2 < s ≤ 3/2, the sum grows as L^{3−2s} and the correction is O(L^{1−2s}) = O(1/L) when s = 1.

The combined rate is $O(L^{1-2\sigma})$ for $1/2 < \sigma \leq 3/2$, improving to $O(1/L^2)$ for $\sigma > 3/2$. ∎

**Remark.** This is a manifestation of Weyl's law [14]: for any self-adjoint operator on a one-dimensional domain of length L with Dirichlet boundary conditions, the spectral zeta function reproduces ζ(2s) through the asymptotic eigenvalue density at the spectral edge. The carry bridge operator is a concrete arithmetic realization of this universal phenomenon.

### 7.3 Numerical verification

Verified numerically for s = 1, 2, 3, 4 (even zeta values) and s = 3/2, 5/2 (odd zeta values ζ(3), ζ(5)):

| s | ζ(2s) | Ratio at L=100 | Ratio at L=500 | Ratio at L=1000 |
|---|-------|----------------|-----------------|-----------------|
| 1 | π²/6 = 1.6449 | 1.6776 | 1.6515 | 1.6482 |
| 2 | π⁴/90 = 1.0823 | 1.1266 | 1.0910 | 1.0867 |
| 3/2 | ζ(3) = 1.2021 | 1.2291 | 1.2075 | 1.2048 |
| 5/2 | ζ(5) = 1.0369 | 1.0542 | 1.0404 | 1.0387 |

All ratios converge to the exact ζ(2s) values at rates consistent with Theorem 3 ($O(L^{1-2\sigma})$ for $\sigma \leq 3/2$, $O(1/L^2)$ for $\sigma > 3/2$).

### 7.4 Universality and non-universality

The ζ(2s) convergence is *universal*: it depends only on the Dirichlet boundary conditions and the length L, not on the specific transition probabilities. Any operator with the same boundary structure produces the same result.

By contrast, the conjectured limit R → −π is *not* universal — it depends on the binary carry structure, the Fejér kernel correlations, and the D-odd constraint (Section 10).

---

## 8. The Bernoulli Hierarchy and the Gaussian Integers

### 8.1 The Bernoulli hierarchy

The carry chain in binary multiplication encodes Bernoulli numbers at successive levels of its spectral decomposition:

| Level | Carry quantity | Value | Bernoulli |
|-------|---------------|-------|-----------|
| 0 | Stationary distribution | 1 | B₀ = 1 |
| 1 | Spectral gap λ₂ | 1/2 | \|B₁\| = 1/2 |
| 2 | Sub-sector ratio α | 1/6 | B₂ = 1/6 |
| 3 | Value cancellation σ₀₀^{c₂=1,c₃=1} | 0 | B₃ = 0 |

**Level 0** is trivial: the carry chain has a unique stationary distribution. **Level 1** is the Diaconis–Fulman theorem [1]: the spectral gap of the carry transfer operator in base b is 1/b, giving 1/2 = |B₁| for binary. **Level 2** is new. The sub-sector ratio α measures the carry population at position D−2, which is directly linked to the boundary layer decomposition of [E, §8.5]: the c_{D−2} = 1 sub-population is precisely the "top" cascade contribution that is blocked in sector (1,0).

**Definition 5 (Sub-sector ratio).** Define α(K) = σ₀₀^{c=1}(K)/σ₀₀^{c=0}(K), where σ₀₀^{c=0} and σ₀₀^{c=1} are the sector-(0,0) partition function restricted to pairs with c_{D−2} = 0 and c_{D−2} = 1, respectively.

**Conjecture 2 (Bernoulli sub-sector ratio).** α(K) → 1/6 = B₂ as K → ∞.

*The evidence is preliminary: the overshoot at K ≥ 20 means the conjecture rests on Richardson extrapolation (3.3 digits) rather than direct convergence. Additional data at higher K would strengthen the case.*

*Numerical evidence.* Exact enumeration gives:

| K | α(K) | $\lvert \alpha - 1/6 \rvert$ | Sig. digits |
|---|------|---------|-------------|
| 7 | 0.15177 | 0.0149 | 1.2 |
| 10 | 0.15958 | 0.0071 | 1.5 |
| 13 | 0.16378 | 0.0029 | 1.9 |
| 16 | 0.16575 | 0.0009 | 2.4 |
| 18 | 0.16652 | 0.0001 | 2.7 |
| 20 | 0.16880 | 0.0021 | 2.0 |
| 21 | 0.16922 | 0.0026 | 2.6 |

The convergence is non-monotonic, with a two-phase structure: approach from below through K = 18, then overshoot at K ≥ 20. The overshoot is present at K = 21 (α = 0.1692, gap = 0.0026). Richardson extrapolation with the ansatz α(K) = α(∞) + a₁K²/2^K + … gives α(∞) ≈ 0.1667 (3.3 digits), consistent with 1/6. The correction structure is K²/2^K rather than pure 1/2^K. The non-monotonic convergence is analogous to the oscillatory behaviour of R_recon in Paper E (§8.3), likely reflecting a Gibbs-type phenomenon in the spectral truncation.

### 8.2 The Euler–Maclaurin connection (speculative)

*This connection is a structural analogy rather than a formal derivation; the numerical coincidence is striking but unproved.*

The carry chain computes ⌊(conv + carry)/2⌋ at each position — a discrete summation with rounding. The Euler–Maclaurin formula for the discrepancy between a discrete sum and its integral is

    Σf − ∫f = B₁·Δf + B₂/2!·Δf' + B₄/4!·Δf''' + ⋯

where B_n are Bernoulli numbers. The carry chain's spectral decomposition mirrors this hierarchy:
- B₁ = −1/2 controls the first-order correction (exponential convergence rate = spectral gap);
- B₂ = 1/6 controls the second-order correction (the sub-sector ratio α, which measures how much of the carry at D − 2 is "rescued" from the c = 1 sub-population).

This connection is structural rather than a formal derivation: the carry chain is a nonlinear discrete dynamical system, not a trapezoidal sum. Nonetheless, the numerical evidence strongly supports the identification.

**Theorem 6 (Cascade Rigidity).** For D-odd sector-(0,0) pairs with $c_2 \in \{0,1\}$:

(i) $\text{carries}[D-1] = \text{carries}[D] = 0$.

(ii) If $c_2 = \text{carries}[D-2] \geq 1$: $M = D-2$ and $\text{val} = c_3 - 1$.

(iii) If $c_2 = 0$ and $c_3 = \text{carries}[D-3] \geq 1$: $M = D-3$ and $\text{val} = c_4 - 1$.

(iv) If $c_2 = c_3 = 0$: $M$ depends on the full carry chain.

*Proof.* In sector (0,0), $x_{K-2} = y_{K-2} = 0$ and $x_{K-1} = y_{K-1} = 1$. The convolution at position $D-2$ is $\text{conv}_{D-2} = x_{K-2} y_{K-1} + x_{K-1} y_{K-2} = 0$. Then $\text{carries}[D-1] = \lfloor(0 + c_2)/2\rfloor = 0$ for $c_2 \in \{0,1\}$, and $\text{carries}[D] = \lfloor(1 + 0)/2\rfloor = 0$. Scanning from $D$ downward, the first nonzero carry is at $D-2$ when $c_2 \geq 1$, giving $M = D-2$ and $\text{val} = \text{carries}[D-3] - 1 = c_3 - 1$. When $c_2 = 0$, the scan continues to $D-3$, etc. □

The theorem extends recursively: at each level $n \geq 2$, if $c_n = \text{carries}[D-n] \geq 1$ while $c_2 = \cdots = c_{n-1} = 0$, then $M = D-n$ and $\text{val} = c_{n+1} - 1$. This recursive extension is verified numerically through $n = 6$ at $K \leq 13$; a general inductive proof for all $n$ is not given.

**Corollary (Level 3: B₃ = 0).** For $c_2 \geq 1$, Theorem 6(ii) gives $\text{val} = c_3 - 1$. When $c_3 = 1$, every pair contributes $\text{val} = 0$, so $\sigma_{00}^{c_2 \geq 1,\, c_3 = 1} = 0$ exactly. This is the carry-chain incarnation of $B_3 = 0$: the corresponding term in the Euler–Maclaurin formula vanishes. Verified for all $K \leq 13$ ($n = 1{,}789{,}942$ pairs at $K = 13$).

The recursive extension produces the same vanishing at every depth: $\sigma_{00}^{c_2 = 0,\, c_3 \geq 1,\, c_4 = 1} = 0$, and so on. At each level, only the "all-zeros" sub-population $(c_2 = c_3 = \cdots = c_n = 0)$ carries non-trivial cascade dynamics.

**Level 4: open.** The non-trivial spectral content telescopes into the $(c_2, c_3, c_4) = (0,0,0)$ and $(0,0,1)$ cells. Seven candidate partition-function ratios, the asymptotic correction structure of $\alpha(K)$, and the spectral gap of the near-boundary carry transfer matrix were tested through $K = 13$; none produces $|B_4| = 1/30$. Identifying the correct observable remains open.

### 8.2a Stopping-time decomposition

Theorem 6, combined with the unit leading carry (Theorem 7, §10.5), yields an exact stopping-time decomposition of the sector ratio:

$$\sigma_{ab} = \sum_{\tau=2}^{D-2} n_{ab}(\tau) \cdot E[\text{val} \mid \tau,\, ab]$$

where $\tau = D - M$ is the depth of the first nonzero carry from the MSB. The Cascade Rigidity Theorem makes this decomposition rigorous: at each $\tau$, val $= c_{\tau+1} - 1$ where $c_{\tau+1} = \text{carries}[D-\tau-1]$.

**Structural exclusion.** In sector $(1,0)$, $\text{conv}_{D-2} = 1$ (Proposition 4), so $\text{carries}[D-2] = 0$ and $P(\tau = 2 \mid \text{sector } 10) = 0$ for all $K$. Sector $(0,0)$ has $\text{conv}_{D-2} = 0$, allowing $\tau = 2$.

**Universality.** Exact enumeration confirms that both $P(\tau \mid ab)$ and $E[\text{val} \mid \tau, ab]$ converge to $K$-independent constants.

**Table 1.** Stopping-time decomposition constants at $K = 21$ (from E45 exact enumeration, $3.4 \times 10^{11}$ D-odd pairs in the $X \leq Y$ half-plane); cross-K stability from $K = 19$ to $K = 21$ confirms convergence to 5–6 significant digits.

| $\tau$ | $P(\tau\mid 00)$ | $P(\tau\mid 10)$ | $E[v\mid\tau,00]$ | $E[v\mid\tau,10]$ | $\Delta(\tau)$ |
|--------|----------------:|----------------:|-----------------:|-----------------:|-------------:|
| 2 | 0.54070 | 0.00000 | −0.00812 | — | — |
| 3 | 0.13713 | 0.45726 | −0.11761 | +0.26971 | +0.38732 |
| 4 | 0.10924 | 0.22875 | −0.04089 | +0.31022 | +0.35111 |
| 5 | 0.08557 | 0.14398 | −0.07338 | +0.33738 | +0.41076 |
| 6 | 0.05007 | 0.07849 | −0.00867 | +0.30657 | +0.31525 |
| 7 | 0.03329 | 0.04308 | −0.01270 | +0.32311 | +0.33580 |
| 8 | 0.01862 | 0.02313 | +0.04186 | +0.33033 | +0.28848 |
| 9 | 0.01169 | 0.01231 | +0.01742 | +0.33307 | +0.31565 |

The sector asymmetry $\Delta(\tau) = E[v\mid\tau,10] - E[v\mid\tau,00]$ is stable to $< 0.00001$ across $K = 19$–$21$ for $\tau \leq 5$, confirming genuine universality. Richardson extrapolation (assuming convergence rate $\rho = 1/2$) yields 7-digit estimates of the limiting constants.

**Identification status.** PSLQ searches against $\{1, \pi, \ln 2, \ln 3\}$ on the Richardson-extrapolated limits produce candidate relations with residuals in the range $10^{-6}$–$10^{-7}$ but no confident closed forms. The precision is sufficient to confirm universality but not to prove algebraic identities. Extending the stopping-time analysis to $K \geq 24$ (computationally feasible with the E45 C implementation) would provide the 10+ digits needed for reliable PSLQ identification.

### 8.2b Failure of the per-position spectral limit

A natural approach to the sector ratio is to analyze the per-position carry transition matrices $T(d)$ directly.  Under D-odd conditioning, the carry chain at each position $d$ induces a stochastic matrix on carry states $\{0, 1, 2, \ldots\}$.  In the Markov (independent-convolution) model, the bulk eigenvalue is $\lambda_2 = 1/2$ at every position; the Linear Mix Hypothesis [E] would predict $\lambda_2 = (1-A^{\ast})/2 + A^{\ast}/2 \approx 0.69$.

Empirical computation of the bulk-averaged transition matrix for $K = 5$–$12$ (experiment P1\_05) reveals a qualitatively different picture:

| K | $\lambda_2(\text{bulk})$ | $A_{\text{eff}}$ |
|---|:---:|:---:|
| 6 | 0.690 | −3.25 |
| 8 | 0.756 | −0.88 |
| 10 | 0.810 | −0.40 |
| 12 | 0.858 | −0.19 |

The second eigenvalue $\lambda_2$ grows monotonically toward 1, and the effective mixing parameter $A_{\text{eff}}$ extracted from the LMH formula diverges to $-\infty$ instead of converging to $A^{\ast} \approx 1.38$.  The per-position transition matrices become *more rigid* (closer to identity) as $K$ increases — the D-odd conditioning progressively constrains the local transitions, which is the dynamical mechanism behind Theorem 6 (Cascade Rigidity).

**Conclusion.**  The Linear Mix Hypothesis cannot hold at the per-position level.  The LMH must describe a *macroscopic* property of the full resolvent

$$(I - \mathcal{K}_{\text{eff}})^{-1}$$

, not a uniform modification of each per-position factor.  This is consistent with the resolvent universality mechanism identified in [E, §8.3]: the spectral resolvent achieves $-\pi$ through collective mode weighting, not through any single position's transition structure.  The stopping-time decomposition (§8.2a) provides the correct analytical framework precisely because it bypasses the per-position transition matrices entirely.

**Boundary-layer concentration.** Approximately 95% of $R(K)$ comes from $\tau \leq 7$. The cumulative ratio through $\tau = 7$ at $K = 14$ is $-2.66$, close to the full $R(14) = -2.62$. The cascade dynamics is overwhelmingly a boundary-layer phenomenon.

**Even-odd alternation.** The successive ratios $P(\tau+1 \mid 00)/P(\tau \mid 00)$ exhibit a clear even-odd oscillation converging toward $1/2$: $\{0.254, 0.796, 0.783, 0.584, 0.664, 0.557, 0.626, \ldots\}$. This is consistent with the Diaconis–Fulman spectral gap $\rho = 1/2$ as the asymptotic rate, modulated by the carry chain's non-stationary structure near the boundary.

### 8.3 The ζ(2) product (preliminary)

*This claim is conditional on two unproved conjectures (α → 1/6 and R → −π); the "new proof of the Basel problem" is speculative until both are established.*

The combination of levels 1 and 2 with the conjectural R = −π gives

    α · R² = (1/6) · π² = ζ(2)

If both α = 1/6 (Conjecture 2) and R = −π (Conjecture 1) can be established analytically from carry arithmetic, this would yield a new proof of the Basel problem ζ(2) = π²/6 — one that begins with integer multiplication rather than Fourier analysis or contour integration.

*Numerical verification.* The product α(K) · R(K)² converges:

| K | α · R² | $\lvert \alpha \cdot R^2 - \zeta(2) \rvert$ |
|---|--------|-----------------|
| 11 | 0.405 | 1.24 |
| 15 | 1.218 | 0.43 |
| 18 | 1.583 | 0.062 |
| 20 | 1.648 | 0.003 |
| 21 | 1.662 | 0.017 |

The product α · R² overshoots ζ(2) at K = 21 because α overshoots 1/6 while R has not yet reached −π. The near-exact agreement at K = 20 (gap = 0.003) is partially accidental: it reflects a cancellation between the α overshoot and the R undershoot at that particular K.

### 8.4 Base specificity

The sub-sector ratio α = 1/6 is **not universal** — it is specific to base 2. Monte Carlo experiments in bases 3, 5, and 7 show no convergence to 1/6 or to any other Bernoulli number. In base 3, α is wildly unstable across K values; in base 5, the sector ratio itself appears rational (R₅ → 5/4) with α₅ far from 1/6. This base specificity makes the binary Bernoulli hierarchy more, not less, remarkable: it reflects a deep structural property of base-2 arithmetic, potentially connected to the ramification of 2 in ℤ[i] (Section 8.6). See [G] for the proof that base 2 is the unique base where the D-parity boundary is a straight line in angular coordinates.

### 8.5 The Euler product perspective

Fulman's multiplicativity property K_a · K_b = K_{ab} for carry transfer operators [3] implies that the Markov spectral sum in base b factors as

    S_b(0) = b²/(b² − 1)

(verified numerically for b = 2, 3, 5, 7, 10). For b = 2 this gives S₂(0) = 4/3. (Note: this spectral zeta sum S₂(0) is a different quantity from the Markov sector ratio R_Markov = +2/3 of Theorem 2; both derive from independent-convolution statistics but measure different observables.) Taking the product over all primes:

    ∏_{p prime} S_p(0) = ∏_p p²/(p² − 1) = ζ(2)

This is a restatement of the classical Euler product, but viewed through carry arithmetic: each prime p contributes its Markov spectral sum factor to ζ(2). The non-Markov (exact) sector ratio introduces L(1, χ₄) = π/4 as the additional factor, connecting to the Dedekind zeta function (Section 8.6).

### 8.6 The Dedekind zeta function of ℚ(i) (speculative)

*This is a heuristic connection: the factorization is classical; attributing both factors to carry arithmetic is conjectural and depends on Conjecture 1.*

The Gaussian integers ℤ[i] have Dedekind zeta function

    ζ_{ℚ(i)}(s) = ζ(s) · L(s, χ₄)

The carry bridge produces:
- ζ(2s) from the spectral zeta function (Theorem 3), via the Markov eigenvalue density
- L(1, χ₄) = π/4 from the sector ratio (Conjecture 1), via the effective χ₄ model (§6.2) and resolvent universality [E]

Both factors of ζ_{ℚ(i)} thus emerge from the same arithmetic object — the carry bridge of binary multiplication. The prime 2 is the unique rational prime that *ramifies* in ℤ[i], satisfying 2 = −i(1+i)². This ramification forces the character χ₄ into the spectral decomposition, explaining at the algebraic level why base 2 is unique [G]: for b ≥ 3, the sector ratio appears to be rational (e.g., R = 5/4 for base 5), and the D-parity boundary is curved in angular coordinates, removing the geometric mechanism that produces L(1, χ₄).

---

## 9. Numerical Evidence

### 9.1 Exact enumeration

The sector ratio R(K) is computed by exhaustive enumeration of all 4^{K−1} digit pairs. Enumeration through K = 21 covers 5.5 × 10¹¹ pairs for K = 21 (C implementation with 10-thread parallelism, ~2.7 hours) [5].

| K | R(K) | $\lvert R(K) + \pi \rvert$ | Sig. digits |
|---|---------|----------|-------------|
| 7 | −0.092 | 3.050 | — |
| 9 | −0.726 | 2.416 | — |
| 11 | −1.554 | 1.588 | — |
| 13 | −2.339 | 0.803 | 0.1 |
| 15 | −2.824 | 0.317 | 0.5 |
| 17 | −3.035 | 0.107 | 1.0 |
| 19 | −3.110 | 0.032 | 1.5 |
| 20 | −3.125 | 0.017 | 1.8 |
| 21 | −3.133 | 0.008 | 2.1 |

All values from the cascade valuation (Definition 1), cross-validated against [E, Table 1].

### 9.2 Richardson extrapolation

The raw convergence R(K) → −π is slow: at K = 21 we have 2.1 significant digits. Richardson extrapolation dramatically improves precision.

**Ansatz 1 (polynomial).** A natural first guess is R(K) = R(∞) + a₁/K² + a₂/K⁴ + ⋯, motivated by the O(1/L²) finite-size correction in the Möbius sum S(A, L) (Proposition 2). However, the true convergence rate is exponential (§9.3), and the polynomial Neville table diverges in practice.

**Ansatz 2 (exponential).** The empirical gap ratios g(K)/g(K−1) converge monotonically to 1/2 (§9.3, Table 1 of [E]), consistent with

    R(K) = R(∞) + c₁ ρ^K + c₂ ρ^{2K} + c₃ ρ^{3K} + ⋯,    ρ = 1/2.

The Neville–Aitken Richardson table with step-size variable h = ρ^K eliminates successive error terms.

**Proposition 3 (Richardson extrapolation).** Applying the Neville–Aitken deferred approach to the limit [9] to the sequence R(K) for K = 7, 9, 11, …, 21 with Ansatz 2 (h = (1/2)^K) yields

    R_extrap = −3.1419… (4.0 significant digits)

using exact R(K) for K ≤ 13 and K ≥ 19, with E44 approximations for K = 14–18. Using only exact values (K ≤ 13, K ≥ 19) the estimate improves to 4.4 significant digits. Leave-one-out cross-validation confirms stability: removing any single K ∈ {7, …, 17} changes the estimate by less than 10⁻⁴.

### 9.3 Convergence rate

The error ε(K) = |R(K) + π| satisfies

    ε(K) = O(K² · 2^{−K})

for K ≥ 11. The consecutive gap ratio g(K)/g(K−1) forms a monotone sequence converging to 1/2 from above: 1.000, 0.954, …, 0.522, 0.491 at K = 21 (see [E, Table 1]). The gap ratio crossing 2 at K = 21 confirms the asymptotic regime. This exponential rate is consistent with the Diaconis–Fulman spectral radius ρ = 1/2 governing the carry chain [4].

*Remark.* A polynomial fit ε ∼ C/K² provides a reasonable approximation for moderate K but misses the asymptotic structure. The exponential characterization ε ∼ C · (1/2)^K captures the true rate. The Richardson extrapolation of §9.2 remains effective because at moderate K the polynomial and exponential forms are numerically similar.

### 9.4 Verification of the closed form

Theorem 1 is verified numerically at 50-digit precision using arbitrary-precision arithmetic. The algebraic identity 4r/(γ(1−r²)) = 2β/(1+2β) holds to full precision for all tested values of β. The formula S(β*, L) → −π is confirmed at L = 1000 with |S + π| < 4 × 10⁻⁵, consistent with the O(1/L²) rate.

---

## 10. The Gap: From S(A) to R(K)

### 10.1 What is proved

Theorem 1 establishes that the abstract alternating sum S(A) = 2(1−A)/(3−2A) is a closed-form Möbius transformation. Setting S = −π yields A* = (3π+2)/(2(π+1)). The spectral zeta function ζ_F(s) → ζ(2s) is proved (Theorem 3). The Dirichlet character χ₄(n) = sin(nπ/2) emerges as the effective spectral weight (§6.2).

### 10.2 The gap and the universality phenomenon

The gap between the abstract sum S(A) and the physical sector ratio R(K) is the central open problem. Three observations constrain its nature:

**Observation 1 (Distinct convergence rates).** The abstract sum converges polynomially: S(A, L) − S(A, ∞) = O(1/L²). The physical sector ratio converges exponentially: R(K) − R(∞) = O(K²/2^K).

**Observation 2 (A_fit does not converge).** Defining A_fit(K) by the equation S(A_fit, 2K−1) = R(K), the sequence A_fit(K) does not converge to A*.

**Observation 3 (Carry projection is exact Markov).** Projecting the exact (non-Markov) transfer operator onto carry degrees of freedom recovers the Markov transfer matrix exactly. The parameter A is therefore *not* a perturbative eigenvalue shift.

These observations rule out deriving A* from the Fejér kernel correlation structure. Nonetheless, both S(A, L) and R(K) converge numerically to the same value −π — a *universality* phenomenon [P2] that we state as a conjecture:

**Conjecture 3 (Universality).** Let $S(A^{\ast}, L) = 2(1-A^{\ast})/(3-2A^{\ast})$ be the abstract Möbius sum (Theorem 1), and let $R(K) = \sigma_{10}(K)/\sigma_{00}(K)$ be the cascade sector ratio (Definition 2). Then

$$\lim_{L \to \infty} S(A^{\ast}, L) \;=\; \lim_{K \to \infty} R(K) \;=\; -\pi,$$

*where the left equality is proved (Theorem 1) and the right is Conjecture 1. The two systems share different microscopic dynamics (polynomial vs exponential finite-size corrections) but the same macroscopic limit.*

The common mechanism is the Dirichlet character χ₄ and the L-function L(1, χ₄) = π/4. [E] identifies **resolvent universality** as the specific mechanism: the spectral resolvent

$$(I - \mathcal{K}_{\mathrm{eff}})^{-1}$$

, applied to the true per-sector cascade profiles, produces a macroscopic sum converging to $-\pi$ through collective mode weighting. The left limit is an exact algebraic evaluation; the right is supported by 4.0-digit exponential Richardson extrapolation through K = 21 (§9).

### 10.3 Structural decomposition

**Proposition 4 (Top-carry constraint).** For a D-odd product n = p·q with D = 2K−1, the convolution at position D−2 equals exactly conv_{D−2} = a + c, where a = p_{K−2} and c = q_{K−2} are the sector bits.

*Proof.* conv_{D−2} = p_{K−2} · q_{K−1} + p_{K−1} · q_{K−2} = a · 1 + 1 · c = a + c. □

**Corollary.** The D-odd condition requires carries[D−1] = 0, hence (a + c + carries[D−2]) < 2. For sector (1,0): a + c = 1, so carries[D−2] = 0, forcing M ≤ D−3. For sector (0,0): carries[D−2] ∈ {0,1}, yielding two sub-populations. This gives σ_{00} = σ_{00}^{top} + σ_{00}^{low}, while σ_{10} = σ_{10}^{low} has no top contribution. [E, §8.5] independently establishes this structural exclusion: P(j=2|sector 10) = 0, confirmed by backward dynamic programming ([E]).

**Proposition 5 (Exact symmetry).** σ_{10}(K) = σ_{01}(K) for all K.

*Proof.* The substitution (p,q) ↦ (q,p) maps sector (1,0) to (0,1), preserving all carries. □

**Theorem 4 (N₁₁ = 0).** For all K ≥ 2, no D-odd pair has sector (1,1).

*Proof.* In sector (1,1), conv_{D−2} = 2. The D-odd condition requires carries[D−1] = 0, but (2 + carries[D−2])/2 ≥ 1, contradiction. □

**Theorem 5 (Closed-form count ratio).** In the continuous limit K → ∞,

$$N_{10}/N_{00} \to \frac{2\ln(4/3) - 1/2}{2\ln(9/8)} = 0.31992784225\ldots$$

with correction O(1/2^K).

*Proof.* Write X = 2^{K−1}(1 + u), Y = 2^{K−1}(1 + v) with u, v ∈ [0, 1). The D-odd condition (1 + u)(1 + v) < 2 defines a region in (u,v). Sectors are determined by ⌊2u⌋ and ⌊2v⌋: area f₀₀ = 2 ln(9/8), f₁₀ = 2 ln(4/3) − 1/2, f₁₁ = 0. Verified to 8 significant digits against exact enumeration (K = 3–21). □

**Corollary (π isolation).** Since R = (⟨val⟩₁₀/⟨val⟩₀₀) × (N₁₀/N₀₀) and N₁₀/N₀₀ involves only ln 2 and ln 3, the entire factor π enters through the carry weight ratio:

$$\langle\text{val}\rangle_{10}/\langle\text{val}\rangle_{00} \to -\pi \cdot \frac{2\ln(9/8)}{2\ln(4/3) - 1/2} = -9.8197\ldots$$

Any proof must explain how the sector bit perturbation modifies the *distribution* of c_{M−1} to produce this ratio. The Dirichlet series approach (character decomposition by n mod 4) does not cleanly factor R, since ~75% of D-odd products are even.

*Remark.* [E, §9.2] falsifies the hypothesis that the *cascade count* ratio n₀₀/n₁₀ converges to π; the limit is algebraic (≈ 25/8). Combined with Theorem 5, this confirms that π is entirely *dynamical*: it resides in the carry weight functional, not in any geometric count.

### 10.4 The carry weight mechanism

The total variation distance TV(j) between the carry distributions in sectors (1,0) and (0,0) reveals a boundary layer with two peaks: TV ≈ 0.16 at j = K−2 (sector bit perturbation) and TV ≈ 0.54 at j = D−2 (top structural constraint), with a plateau TV ≈ 0.06 in between.

Four steps link the sector bit to the carry weight val = c_{M−1} − 1:

1. **Convolution perturbation (position K−2).** The sector bit a = 1 adds q_{j−K+2} to conv_j for j ∈ [K−2, 2K−3], generating additional carries in the bulk.

2. **Carry propagation.** The extra carries propagate upward through ~K positions; the mean carry at each position is ~0.35 higher in sector (1,0) [F].

3. **Top constraint (position D−2).** The D-odd condition forces carries[D−2] = 0 in sector (1,0), eliminating the "top" population and shifting M to D−3 or lower.

4. **Val distribution shift.** The higher carry levels in sector (1,0), combined with the constrained top, produce 30% with c_{M−1} = 2 (val = +1) vs 5.5% in sector (0,0)_{low}, driving σ_{10} > 0 while σ_{00} < 0.

**Spectral connection.** The sector bit perturbation is a ramp (constant ⟨Δconv_j⟩ = 1/2) over half the chain. In the bridge eigenmode basis sin(nπj/L), the ramp Fourier coefficients Δv_n = [cos(nπ/2) − cos(nπ)]/(nπ) combined with the source-point factor sin(nπ/2) satisfy:

$$\sin(n\pi/2) \cdot [\cos(n\pi/2) - \cos(n\pi)]
= \sin(n\pi/2)\ \text{if } n \text{ is odd},\ \text{and } 0\ \text{if } n \text{ is even}.$$

The even-mode contributions cancel exactly; the odd-mode sum becomes Σ sin(nπ/2)/n = π/4 = L(1, χ₄). This explains conceptually how the Leibniz series enters the problem.

**Remaining normalization gap.** With the Diaconis–Fulman spectral gap λ = 1/2 taken as constant, the spectral sum evaluates to 1/2, not −π. The discrepancy factor −2π must originate from: (i) the mode-dependent eigenvalue spectrum of the non-Markovian carry chain; (ii) the non-local observation functional val = c_{M−1} − 1 evaluated at the random stopping position M−1; and (iii) the redistribution of the M = D−2 population (~55% of sector (0,0)) into M ≤ D−3 in sector (1,0). All three concern carry weights, consistent with the π-isolation corollary. The analytical resolution is identified in [E, §9] as the **Holte bridge / Doob h-transform** matching between the bulk carry distribution (conditioned on D-odd) and the boundary layer exclusions.

---

### 10.5 Analytical evidence: continuum limit and exact symbolic forms

**Carry recursion from MSB.** Define carry$[d]$ = carry at position $D-1-d$ from below, with recursion carry$[d] = 2 \cdot$carry$[d-1] + b_d - C_d$ and carry$[0] = 0$.

**Theorem 7 (Unit leading carry).** Let $d_0$ be the smallest $d \geq 1$ with carry$[d] \geq 1$. Then carry$[d_0] = 1$ exactly.

*Proof.* Since carry$[d_0-1] = 0$, carry$[d_0] = b_{d_0} - C_{d_0}$. For carry $\geq 1$ with $b_{d_0} \in \{0,1\}$ and $C_{d_0} \geq 0$: necessarily $C_{d_0} = 0$ and $b_{d_0} = 1$. □

**Corollary (Universal val formula).** $\text{val} = 1 + b_{d_0+1} - C_{d_0+1}$ for all depths $d_0$ and both sectors.

The continuum function $\text{val}_\infty(u,v) = \lim_{K \to \infty} \text{val}(u,v;K)$ is a well-defined step function taking values in $\{-1,0,+1\}$, giving sector sums $\sigma_{ac} = \int\!\!\int_{\text{sector } ac,\, \text{D-odd}} \text{val}_\infty(u,v)\, du\, dv$.

**Exact symbolic forms.** Each $\sigma(d_0)$ admits an exact closed form as a rational linear combination of $\{1, \ln p : p \text{ prime}, p \leq 2^{d_0+2}\}$:

$$\sigma_{00}(1) = \tfrac{1}{2} + \tfrac{9}{2}\ln 2 - 2\ln 3 - 3\ln 5 + \tfrac{7}{4}\ln 7$$

The individual coefficients grow as $\sim 4^{d_0}$ while $\sigma(d_0) = O(2^{-d_0})$. After 7 depth levels, the cumulative $\ln 2$ coefficient reaches ~1523, while the total $\sigma_{00} \approx -0.007$ — catastrophic cancellation (ratio $> 10^5$) essential to the emergence of $\pi$. At each depth $d_0$, new primes in $[2^{d_0+1}, 2^{d_0+2}]$ enter, so the infinite series involves **all primes**. By Baker's theorem, no finite combination of logarithms of algebraic numbers equals $\pi$; the infinite series circumvents this obstruction.

| $d_{\max}$ | $\sigma_{10}/\sigma_{00}$ | $\lvert R + \pi \rvert$ | digits |
|-----|-----|-----|-----|
| 5 | $-2.695$ | $4.5 \times 10^{-1}$ | 0.85 |
| 10 | $-3.098$ | $4.3 \times 10^{-2}$ | 1.86 |
| 12 | $-3.129$ | $1.2 \times 10^{-2}$ | 2.41 |
| 14 | $-3.139$ | $2.1 \times 10^{-3}$ | 3.18 |
| 15 | $-3.14146$ | $1.4 \times 10^{-4}$ | **4.36** |

Beyond $d_0 = 15$, catastrophic cancellation exceeds float64 precision. Extended precision (60 digits; [E]) reaches $d_0 = 18$ but reveals a structural limitation: the depth-first carry-free enumeration captures only 17–33% of the D-odd configuration space, causing the ratio to worsen beyond $d_0 = 15$. This is a method limitation, not a failure of the $-\pi$ conjecture.

## 11. Open Problems

1. Prove $\sigma_{10}/\sigma_{00} = -\pi$ analytically. The decomposition $R = R_0 + \Delta R$ (where $R_0 = -3.9312\ldots$ involves only $\ln 2$ and $\ln 3$, and $\Delta R = +0.7896\ldots$ carries the entire transcendental content) reduces the problem to showing that the cascade correction $\Delta R = |R_0| - \pi$. Six equivalent formulations are catalogued in §10.5; the depth series $\sum_{d=1}^\infty \sigma(d)$ converges to 4.36 significant digits at $d_{\max} = 15$. The most promising analytical path is the **stopping-time series** (§8.2a): derive $P(\tau \mid ab)$ and $E[\text{val} \mid \tau, ab]$ analytically from the carry chain structure, then sum $R = \sum_\tau [n_{10}(\tau) \cdot E[\text{val} \mid \tau, 10] + \cdots]$. The sector asymmetry $\Delta(\tau) = E[v \mid \tau, 10] - E[v \mid \tau, 00]$, which converges to $K$-independent constants (Table 1, §8.2a), is the core quantity: it encodes the entire transcendental content of R. The per-position spectral approach (Doob h-transform of local transition matrices) is ruled out: §8.2b shows that the per-position eigenvalues grow with $K$ instead of stabilizing, so the LMH must operate at the global resolvent level, not locally. Current status from the Step3 bridge: the `χ₄ -> L(1,χ₄)` part is reduced to a scalar closure `R(∞)=C·L(1,χ₄)` with data-constrained `C∈[-4.0032,-3.9896]`, and the dominant decay rate is reduced to proving the analytic `1/2` bound for the conditioned chain. See [E, §9] for the updated closure map.

2. **Base universality.** For base b ≥ 3, the Markov sum S_b(0) = b²/(b²−1). Is the exact sector ratio R_b rational for all b ≥ 3? Exact enumeration for base 3 and base 5 would test this. The conjecture R₅ = 5/4 has preliminary numerical support.

3. **Higher Bernoulli levels.** The hierarchy B₀, |B₁|, B₂ is confirmed through Level 2. Level 3 (B₃ = 0) is proved by Theorem 6: the cascade rigidity forces $\sigma_{00}^{c_2 \geq 1,\, c_3=1} = 0$ exactly, and the recursive extension produces the same vanishing at every depth. The Level 4 signal |B₄| = 1/30 does not appear as a simple partition-function ratio, asymptotic correction coefficient, or spectral gap (all tested through K = 13). The non-trivial dynamics telescopes into the $(c_2, c_3, \ldots) = (0, 0, \ldots)$ sub-population at each level. Identifying the correct Level 4 observable remains open.

4. **Catalan's constant from carries.** The value L(2, χ₄) = G = 0.9159… (Catalan's constant) is the next L-function value after L(1, χ₄) = π/4. Can it be extracted from carry chain data?

5. **ζ(s) at odd arguments** (speculative). The Dedekind factorization ζ(s) = ζ_{ℚ(i)}(s)/L(s, χ₄) gives ζ(s) at all arguments if L(s, χ₄) is accessible from carries for general s. The case ζ(3) = 1.20206… (Apéry's constant) is the most important test.

6. **New proof of the Basel problem** (conditional). If both α = 1/6 (Conjecture 2) and R = −π (Conjecture 1) are proved from carry arithmetic, then ζ(2) = αR² = π²/6 follows. This depends on two currently unproved conjectures; neither is close to resolution.

---

## References

1. P. Diaconis and J. Fulman, *Carries, shuffling, and symmetric functions*, Adv. Appl. Math. **43** (2009), 176–196.
2. P. Diaconis and J. Fulman, *Carries, shuffling, and an amazing matrix*, Amer. Math. Monthly **116** (2009), 788–803.
3. J. Fulman, *The carries process revisited*, preprint (2023).
4. J. M. Holte, *Carries, combinatorics, and an amazing matrix*, Amer. Math. Monthly **104** (1997), 138–149.
5. D. E. Knuth, *The Art of Computer Programming, Volume 2: Seminumerical Algorithms*, 3rd ed., Addison-Wesley, 1997. (Chapter 4.3: Multiple-precision arithmetic.)
6. L. Euler, *De summis serierum reciprocarum*, Comment. Acad. Sci. Petrop. **7** (1740), 123–134.
7. R. Remmert, *Classical Topics in Complex Function Theory*, Graduate Texts in Mathematics **172**, Springer, 1998. (Chapter 9: The number π.)
8. M. V. Berry and J. P. Keating, *The Riemann zeros and eigenvalue asymptotics*, SIAM Rev. **41** (1999), 236–266.
9. L. F. Richardson and J. A. Gaunt, *The deferred approach to the limit*, Phil. Trans. R. Soc. A **226** (1927), 299–361.
10. I. S. Gradshteyn and I. M. Ryzhik, *Table of Integrals, Series, and Products*, 7th ed., Academic Press, 2007. (Entry 1.421.3.)
11. P. G. L. Dirichlet, *Beweis des Satzes, dass jede unbegrenzte arithmetische Progression…*, Abh. Königl. Preuss. Akad. Wiss. (1837), 45–81.
12. G. Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*, 3rd ed., Graduate Studies in Mathematics **163**, AMS, 2015. (Chapter II.5: Selberg–Delange method.)
13. T. Apostol, *Introduction to Analytic Number Theory*, Springer, 1976. (Chapters 6, 12: Dirichlet characters and L-functions.)
14. H. Weyl, *Über die asymptotische Verteilung der Eigenwerte*, Nachr. Ges. Wiss. Göttingen (1911), 110–117.
15. D. Zagier, *Values of zeta functions and their applications*, in *First European Congress of Mathematics*, Progress in Math. **120**, Birkhäuser, 1994, 497–512.
16. E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed. (revised by D. R. Heath-Brown), Oxford, 1986.
17. [P2] Companion paper: *The sector ratio in binary multiplication: from Markov failure to transcendence*.
18. [G] Companion paper: *The angular uniqueness of base 2 in positional multiplication*.
19. [E] Companion paper: *The Trace Anomaly of Binary Multiplication*. (Identifies conditionally the mechanism for $R \to -\pi$ via resolvent universality, assuming the LMH.)
20. [A] Companion paper: *Spectral theory of carries in positional multiplication*. (Foundation: the $m$-bit Equidistribution Lemma extending Diaconis–Fulman to the transfer operator.)
21. [F] Companion paper: *Exact covariance structure of binary carry chains*. (Carry expectation $E[c_j] = (j-1)/4$, off-diagonal covariance $\text{Cov}(c_j, g_i h_{j-i}) = 1/8$.)
22. [B] Companion paper: "Carry Polynomials and the Euler Product: An Approximation Framework," this series.
