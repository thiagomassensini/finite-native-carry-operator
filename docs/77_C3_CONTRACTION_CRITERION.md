# C3 contraction ladder and exact-zero criterion

## What the finite ladder certifies

The oriented certificates at cutoffs \(M=16\,364\), \(M=65\,536\), and
\(M=131\,072\) enclose the same limiting C3 operator.  Each refined interval
is contained in the preceding interval, and every certificate proves
existence and uniqueness of a stationary point for the same limiting-energy
equation.  Therefore all three certificates bound one and the same point
\(t_*\).

The machine-checkable ladder is
[`certificates/c3/c3_oriented_contraction_ladder.json`](../certificates/c3/c3_oriented_contraction_ladder.json).
Its certified bounds are:

| Cutoff | Interval radius | \(\lVert R_\infty(t_*)\rVert\) upper bound | Energy upper bound |
|---:|---:|---:|---:|
| \(16\,364\) | \(5\times10^{-16}\) | \(2.510236\times10^{-15}\) | \(6.301281\times10^{-30}\) |
| \(65\,536\) | \(2\times10^{-20}\) | \(1.661316\times10^{-19}\) | \(2.759971\times10^{-38}\) |
| \(131\,072\) | \(4\times10^{-22}\) | \(3.477616\times10^{-21}\) | \(1.209382\times10^{-41}\) |

The norm bound contracts by a factor greater than \(1.5109\times10^4\), and
the corresponding energy bound contracts by its square, greater than
\(2.2830\times10^8\).  The oriented resultant-tail remainder contracts by a
factor greater than \(2.0612\times10^3\).  Every energy upper bound in the
ledger is stored as the exact integer square of the corresponding norm upper
bound, including its decimal exponent.

These are three rigorously verified members of a finite ladder.  Their strong
contraction is evidence for the next construction, but finitely many members
alone do not prove that an infinite family converges to zero.

## Exact logical bridge

For one fixed limiting point \(t_*\), suppose a nonnegative family \(e_M\)
satisfies

\[
\lVert R_\infty(t_*)\rVert^2\le e_M
\quad\text{for every }M,
\qquad e_M\longrightarrow0.
\]

Then the fixed nonnegative quantity on the left is zero, hence
\(R_\infty(t_*)=(0,0)\).  This implication is formalized without additional
assumptions by
`VanishingLimitResidualCertificate.witness_zero` in
[`FiniteNativeCarryOperator/Certification/Contract.lean`](../FiniteNativeCarryOperator/Certification/Contract.lean).

## Remaining certification obligation

The missing mathematical step is now precise: construct a cutoff-uniform
family of oriented interval certificates for the same \(t_*\), derive an
explicit energy upper bound \(e_M\), and prove \(e_M\to0\).  More isolated
cutoffs can test and guide that formula, but they do not replace the uniform
proof.

Accordingly, the ladder records both
`infinite_vanishing_bound_family_certified: false` and
`limiting_vector_zero_certified: false`.

## Reverification

```bash
python3 scripts/verify_c3_contraction_ladder.py \
  certificates/c3/c3_oriented_contraction_ladder.json --recompute
```

The verifier rebuilds the ladder from the three oriented certificates, checks
their hashes, proves interval nesting, checks every strict contraction, and
validates the exact square relation between norm and energy bounds.
