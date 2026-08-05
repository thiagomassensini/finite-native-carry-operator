# Finite Green, boundary, and return audit

The scanner contains an exact finite reconstruction check for an arbitrary
real-plane state sequence \(f\):

\[
f = G(Bf)+R(\operatorname{Tr}f).
\]

Here `B` is the second-difference curvature, `G` is double summation with zero
initial trace, and the return term is the affine sequence carrying the initial
position and slope.  The numerical audit reports reconstruction, curvature,
trace, and polarized energy-ledger errors.

This finite identity is independent of the asymptotic corrector.  The contextual
Markov note is retained with its original firewall: emission/return on a divisor
tree does not automatically become a Green identity on an edge space without
an explicit intertwiner or boundary synthesis.
