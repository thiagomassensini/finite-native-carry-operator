# Dyadic asymptotic tail corrector

## Native tail scale

For the native real state, the centered second difference starts with a second
derivative term.  Summing its tail produces a leading real rotational layer

\[
M^{-3/2}\operatorname{Rot}(-t\log M)v_0(t).
\]

When \(M\) doubles, this term contracts by \(2^{-3/2}\) and rotates by phase
\(-t\log2\).  Subsequent layers have radial scales
\(M^{-5/2}, M^{-7/2},\ldots\), each carrying the same real-plane rotation law.

## Annihilator

For one layer with contraction \(\rho\) and phase \(\theta\), the real sequence
is annihilated by

\[
L^2-2\rho\cos(\theta)L+\rho^2.
\]

The canonical corrector multiplies the factors for \(\rho_0=2^{-3/2}\) and
\(\rho_1=2^{-5/2}\), giving a fourth-order recurrence.  Phase consistency is
solved by fixed-point iteration, not Newton.

## Holdout discipline

Each rolling holdout predicts a cutoff that was excluded from its training
window.  The final preserved errors for \(M=2^{21}\) and \(M=2^{22}\) are
approximately \(1.47\times10^{-16}\) and \(5.18\times10^{-17}\), respectively.

The consensus center and its radius remain asymptotic numerical evidence.  A
rigorous result requires an explicit bound on all omitted linear and nonlinear
tail terms.
