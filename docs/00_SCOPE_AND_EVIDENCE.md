# Scope and evidence contract

## Object under study

The repository studies the finite real-plane operator obtained by applying
explicit centered-bracket cameras to

\[
\psi_t(n)=n^{-1/2}(\cos(-t\log n),\sin(-t\log n)).
\]

The finite cutoff \(M\), the camera geometry, the seed set, and every bracket
radius are part of the object and are always recorded.

## Evidence levels

| Mark | Meaning |
|---|---|
| **[D]** | definition implemented identically in Python and represented in Lean |
| **[F]** | exact finite identity or theorem |
| **[N]** | reproducible numerical observation at named precision and cutoff |
| **[A]** | asymptotic model validated by holdouts but not rigorously remainder-bounded |
| **[C]** | interval-certified statement with machine-checkable enclosure and tail proof |

The reconstructed working tree contains [D], [F], [N], and [A] material plus
several linked [C] layers: the fixed-cutoff C3 minimum at \(M=16\,364\),
uniform C3 derivative tails, a unique strict minimum of the limiting energy,
and three nested oriented refinements.  No limiting vector zero has yet been
certified.  At that minimum the remaining vector-zero question has been
reduced exactly to the vanishing of the oriented real-plane determinant
against the certified nonzero velocity.

An oriented boundary-jet refinement further encloses that limiting minimum in
a radius \(5\times10^{-16}\) interval and places a rigorous
\(2.510236\times10^{-15}\) upper bound on the resultant norm at the stationary
point.  The positive upper bound is not recorded as an exact zero.

Two further oriented certificates at \(M=65\,536\) and \(M=131\,072\) are
nested inside the first and contract that bound to
\(1.661316\times10^{-19}\) and \(3.477616\times10^{-21}\), respectively.
Their finite contraction ladder proves that all three certificates enclose
the same limiting stationary point and records exact squared energy bounds.

The residual decomposition proves uniformly that the analytic tail term tends
to zero like at least \((M+1)^{-5}\).  It leaves two explicit obligations:
prove that the corrected-center residual \(Q_M\to0\), and construct nested
radii with \(r_M\to0\).  The implication from these three component limits
to an exact real-plane zero is formalized in Lean.

The stationary-localization refinement derives every finite radius from the
corrected stationary residual and one positive lower bound for
\(H_\infty'\).  It constructs the unique corrected-root family from the
certified cutoff threshold and proves that its derived radius tends to zero.

## Nonclaims

- A local minimum of \(\lVert R_M(t)\rVert^2\) is not silently renamed a vector
  zero of \(R_M\).
- Agreement among multiprecision stages at fixed \(M\) does not imply agreement
  with \(M\to\infty\).
- The empirical radius of the asymptotic corrector is not an interval enclosure.
- Finitely many contracting certified cutoffs do not by themselves prove convergence of
  the certified upper bounds to zero.
