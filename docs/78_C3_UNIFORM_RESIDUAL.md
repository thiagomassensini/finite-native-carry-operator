# C3 cutoff-uniform residual decomposition

## Purpose

The finite contraction ladder bounds the limiting resultant at one shared
stationary point, but a decreasing list of positive bounds is not yet an
infinite convergence proof.  This document separates each bound into the
pieces that do and do not already have cutoff-uniform control.

Let \(t_*\) be the unique limiting stationary point, and let
\(I_M=[c_M-r_M,c_M+r_M]\) be a certified interval containing it.  Write

\[
R_\infty(t)=R_M(t)+J_M(t)+E_M(t),
\]

where \(J_M\) is the fifth-order oriented boundary jet and
\(\lVert E_M\rVert\le\eta_M\).  Define

\[
Q_M=\lVert R_M(c_M)+J_M(c_M)\rVert.
\]

If one constant \(V\) bounds \(\lVert R_\infty'(t)\rVert\) on the largest
certified interval, then interval nesting and the mean-value bound give

\[
\boxed{
\lVert R_\infty(t_*)\rVert
\le Q_M+\eta_M+Vr_M.
}
\]

The ledger
[`c3_uniform_residual_decomposition.json`](../certificates/c3/c3_uniform_residual_decomposition.json)
stores all three terms independently.  It uses the single certified cap

\[
V<4.628890
\]

from the largest interval, so the localization term is uniform across every
nested entry.

## Analytic tail component

For a common bound \(|t|\le T\), set

\[
P_6(T)=\prod_{j=0}^{5}\sqrt{T^2+(j+\tfrac12)^2}.
\]

The oriented resultant remainder has the explicit bound

\[
\eta_M(T)=P_6(T)\left[
\frac7{660}\bigl(3(M+1)\bigr)^{-11/2}
+\frac1{5940}(3M-1)^{-11/2}
\right].
\]

For every integer \(M\ge2\), both lower arguments are at least \(M+1\).
Since \(x^{-11/2}\le x^{-5}\) for \(x\ge1\),

\[
\eta_M(T)
\le
P_6(T)\frac{16}{1485}(M+1)^{-5}.
\]

Thus \(\eta_M(T)\to0\) for fixed \(T\).  This is a genuine
cutoff-uniform conclusion, not an inference from the three computed entries.
The polynomial witness is intentionally simple and looser than the sharp
remainder.  Lean formalizes both the convergence of the polynomial envelope
and the squeeze step for any nonnegative error below it in
`tendsto_zero_of_le_polynomialTailEnvelope`.

## Certified finite decomposition

| \(M\) | \(Q_M\) | \(\eta_M\) | \(Vr_M\) | Total norm envelope |
|---:|---:|---:|---:|---:|
| \(16\,364\) | \(4.43601\times10^{-17}\) | \(1.07184\times10^{-16}\) | \(2.31445\times10^{-15}\) | \(2.46599\times10^{-15}\) |
| \(65\,536\) | \(1.76687\times10^{-23}\) | \(5.19980\times10^{-20}\) | \(9.25778\times10^{-20}\) | \(1.44594\times10^{-19}\) |
| \(131\,072\) | \(1.06214\times10^{-24}\) | \(1.14906\times10^{-21}\) | \(1.85156\times10^{-21}\) | \(3.00167\times10^{-21}\) |

The third interval has radius \(4\times10^{-22}\).  Its core residual
contracts by a factor greater than \(16.63\) relative to the preceding
entry.  Across all three entries, the core residual and radius both contract
strictly.

## Exact remaining obligations

The decomposition identifies two remaining infinite statements:

1. construct corrected centers \(c_M\) uniformly and prove \(Q_M\to0\);
2. construct the nested radii uniformly and prove \(r_M\to0\).

The tail component is complete.  If the two remaining statements are supplied,
then \(Q_M+\eta_M+Vr_M\to0\), and exact vector vanishing follows.

The radius implication is refined in
[`79_C3_STATIONARY_LOCALIZATION.md`](79_C3_STATIONARY_LOCALIZATION.md):
a common positive stationary-slope margin produces an explicit vanishing
radius for any cutoff-uniform family of exact corrected stationary roots.

This three-component implication is formalized in Lean as
`DecomposedVanishingLimitResidualCertificate.witness_zero`.  The present
ledger correctly retains
`core_residual_Q_M_infinite_vanishing_family_certified: false`,
`radius_infinite_vanishing_family_certified: false`, and
`limiting_vector_zero_certified: false`.

## Reverification

```bash
python3 scripts/verify_c3_uniform_residual.py \
  certificates/c3/c3_uniform_residual_decomposition.json --recompute
```

The verifier rebuilds the finite sums, boundary jets, sharp remainders,
polynomial witnesses, common velocity cap, localization terms, and energy
envelopes using directed rounding.
