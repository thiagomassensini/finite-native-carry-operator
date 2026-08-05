# Roadmap: real interval calculation and certification

The certification line replaces empirical localization claims by explicit
machine-checkable real enclosures.  The path is deliberately staged.

## Phase A — directed-rounding finite evaluator

1. Implement outward-rounded enclosures for `log`, `sin`, `cos`, powers, and
   every real-plane addition.
2. Enclose each seed, bracket coordinate, resultant coordinate, emitted energy,
   and derivative on a closed interval in `t`.
3. Cross-check two independent real backends, for example Arb/FLINT and MPFI
   or MPFR-based Rust.
4. Serialize every enclosure with precision, rounding mode, source hash, and
   exact camera/cutoff metadata.

Status: implemented with Arb real balls for the first C3 certificate.  The
second independent backend remains open.

## Phase B — finite stationary minimum

1. Enclose \(R_M\cdot R_M'\) on subintervals and prove a unique sign change.
2. Enclose the second derivative of \(\lVert R_M\rVert^2\) away from zero.
3. Produce a `FiniteMinimumCertificate` matching the Lean contract.
4. Separately test whether either coordinate enclosure contains zero; do not
   identify a stationary minimum with a vector zero.

Status: completed for C3 at \(M=16\,364\) on the certified interval.  The
certificate proves a unique strict finite minimum and excludes a finite vector
zero there.

## Phase C — finite vector zero, when applicable

A one-parameter curve in \(\mathbb R^2\) needs an internal, problem-specific
argument.
Acceptable routes include a proven common scalar factor, an exact phase
constraint reducing the system to one scalar equation, or a validated
transversality/topological argument.  Merely observing a tiny norm is not enough.

No external comparison theory is admissible in this phase.

Status: for C3, the limiting velocity is certified nonzero throughout the
limiting-minimum interval.  The real-plane Lagrange identity therefore reduces
the vector-zero question at the unique stationary point to the single scalar
condition \(\det(R_\infty,R'_\infty)=0\).  Certification of that determinant
condition remains open.

## Phase D — rigorous real tail

1. Derive Euler–Maclaurin or summation-by-parts enclosures for every camera
   channel.
2. Bound the remainder uniformly on the certification interval.
3. Prove that the finite-to-limit perturbation is smaller than the uniqueness
   margin.
4. Transport the finite certificate to a `LimitZeroCertificate` or to a
   certified limiting minimum, according to what the enclosed equations support.

Status: completed for the C3 resultant and its first two time derivatives on
the certified domain at \(M=16\,364\).  The positive-kernel tail certificate,
combined with a ten-cell Arb cover, proves a unique strict minimum of the
limiting C3 energy.  It does not yet prove a limiting vector zero.

Refinement status: the oriented fifth-order boundary jet is also certified.
Its sixth-derivative remainder relocalizes the unique limiting minimum to a
radius \(5\times10^{-16}\) interval and bounds the resultant norm there by
\(2.510236\times10^{-15}\).  Equality to zero remains a separate obligation.
The independent \(M=65\,536\) reinforcement contracts these figures to radius
\(2\times10^{-20}\) and norm bound \(1.661316\times10^{-19}\).
The \(M=131\,072\) reinforcement contracts them again to radius
\(4\times10^{-22}\) and norm bound \(3.477616\times10^{-21}\).

The resulting finite contraction ladder verifies that all three nested
intervals enclose the same limiting stationary point.  The residual
decomposition proves the analytic tail component tends to zero with an
explicit \((M+1)^{-5}\) witness and uses one velocity cap for every nested
interval.

The stationary-localization ledger further derives
\(r_M=(|h_M(c_M)|+\delta_M)/m\) from a common positive slope margin.  For
ideal corrected roots \(h_M(c_M)=0\), both \(\delta_M\) and \(r_M\) have
explicit polynomial limits zero.  The ledger constructs that root family from
its certified threshold; the remaining Phase D obligation is to prove
\(Q_M\to0\) along it.  The finite ladder alone is not an exact-zero proof.

## Phase E — Lean ingestion

- Parse the certificate ledger into exact rationals or dyadics.
- Prove enclosure composition lemmas once.
- Verify each generated witness without `axiom` or `sorry`.
- Keep source digests and release commit in the theorem registry.

Status: the abstract final implication is complete.  In addition to
`VanishingLimitResidualCertificate`, Lean now contains
`DecomposedVanishingLimitResidualCertificate`.  Its theorem `witness_zero`
proves exact real-plane vanishing when the corrected-center, oriented-tail,
and localization components all tend to zero.  Ingestion of the two remaining
uniform C3 component proofs into that contract remains open.

The additional Lean structure `StationaryLocalizationCertificate` turns a
vanishing stationary-equation error and a positive slope margin into
vanishing witness-to-center radii.  This completes the abstract localization
bridge used by the new ledger.
