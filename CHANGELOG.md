# Changelog

## Unreleased

- Reconstruct the working repository root from the audited finite-operator
  sources while preserving the raw staging directory locally.
- Add a real Arb directed-rounding evaluator and exact-integer certificate
  verifier.
- Certify a unique strict C3 finite minimum at `M = 16,364` on a decimal-radius
  `1e-40` interval and exclude a finite vector zero from that interval.
- Add explicit positive-kernel C3 tail bounds for the resultant and its first
  two time derivatives, and certify a unique strict minimum of the limiting C3
  energy on a decimal-radius `5e-5` interval.
- Certify nonvanishing limiting velocity on that interval and reduce the C3
  vector-zero question at its stationary point to one oriented real
  determinant, with the equivalence formalized in Lean.
- Add a fifth-order oriented C3 boundary jet with a sixth-derivative interval
  remainder, relocalize the limiting minimum to radius `5e-16`, and bound its
  resultant norm by `2.510236e-15` without claiming exact vanishing.
- Reinforce the same oriented limit certificate at `M = 65,536`, contracting
  the stationary interval radius to `2e-20` and its resultant bound to
  `1.661316e-19`.
- Add a third oriented reinforcement at `M = 131,072`, with interval radius
  `4e-22` and resultant bound `3.477616e-21`.
- Add a machine-checkable nested contraction ladder for the three oriented C3
  certificates, including exact squared energy bounds and explicit finite-only
  nonclaims.
- Decompose the stationary residual as `Q_M + eta_M + V*r_M`; prove the
  analytic tail component has an explicit uniform `(M+1)^-5` witness, and
  isolate the two remaining infinite obligations.
- Add a stationary-localization ledger using one positive limiting slope
  margin, certify derived radii inside all three oriented intervals, and prove
  a polynomial vanishing-radius witness for ideal corrected stationary roots.
- Formalize the cutoff-uniform vanishing-residual criterion in Lean: energy
  upper bounds tending to zero at one fixed limiting witness force exact
  real-plane vanishing.
- Formalize the three-component residual criterion in Lean so the
  corrected-center, tail, and localization limits remain separate.
- Formalize stationary localization by error divided by a positive slope,
  including witness containment and convergence of the velocity-times-radius
  term.
- Express the asymptotic corrector entirely through real-plane contraction and
  rotation.
- Confirm structurally that C2 is the radius-one sector of C4 and numerically
  that cross-camera finite-minimum dispersion contracts with cutoff.
- Repair the finite-energy faithfulness proof so the Lean root builds with
  warnings treated as errors.

## v0.1.0 — 2026-08-05

- Establish the finite native carry operator as an independent real-plane
  laboratory object.
- Preserve C2 and all natural saturated camera geometries, including even
  antipodal channels.
- Add CPU/CUDA discovery and mechanism diagnostics.
- Add arbitrary-precision Ridder localization without Newton.
- Archive C3 cutoff ledgers through `M = 4,194,304`.
- Add one-, two-, and three-layer dyadic tail reports with rolling holdouts.
- Add Lean finite definitions and interval-certification contracts.
- Add citation, Zenodo, provenance, checksums, CI, and an immutable release gate.
