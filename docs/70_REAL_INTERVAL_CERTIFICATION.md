# Real interval certification

## Scope

The certificate is built entirely in the native real plane.  Write

\[
R_M(t)=(X_M(t),Y_M(t)),\qquad
F_M(t)=X_M(t)^2+Y_M(t)^2.
\]

The scalar stationary numerator and its derivative are

\[
H_M=X_MX'_M+Y_MY'_M=\frac12F'_M,
\]

\[
H'_M=(X'_M)^2+(Y'_M)^2+X_MX''_M+Y_MY''_M.
\]

Only real directed-rounding balls are used.  No external comparison theory or
non-real coordinate participates in the construction or its interpretation.

## Certified finite result

For camera 3 and cutoff \(M=16\,364\), the ledger
[`certificates/c3/c3_m16364_raw_minimum.json`](../certificates/c3/c3_m16364_raw_minimum.json)
uses the interval centered at the recorded multiprecision minimum with decimal
radius \(10^{-40}\).  It proves:

1. \(H_M\) is strictly negative at the left endpoint;
2. \(H_M\) is strictly positive at the right endpoint;
3. \(H'_M\) is strictly positive on the full interval;
4. both coordinate enclosures of \(R_M\) exclude zero on the full interval.

Continuity and the first three facts give one and only one stationary point in
the interval, and it is a strict finite minimum of \(F_M\).  The fourth fact
proves separately that the interval contains no finite vector zero.

## Machine-checkable representation

Every ball is serialized as three integers:

```text
(midpoint_integer +/- radius_integer) * 10^exponent10
```

The verifier checks signs using exact integer arithmetic, recomputes every
source SHA-256 digest, and confirms that the published claims follow from the
enclosures:

```bash
python3 scripts/verify_interval_certificate.py \
  certificates/c3/c3_m16364_raw_minimum.json --recompute
```

With `--recompute`, the verifier repeats every directed-rounding evaluation and
then compares the rebuilt enclosures, conditions, claims, versions, and source
digests with the stored ledger.

## Exact logical boundary

This certificate establishes a fixed-cutoff minimum and excludes a fixed-cutoff
vector zero in its interval.  It does not certify a limiting zero.  That next
claim requires uniform real tail bounds for the resultant and its derivatives,
followed by a separate existence argument for simultaneous coordinate
vanishing.

The uniform C3 tail bounds and the transport to a unique limiting stationary
minimum are now established in
[`75_C3_REAL_TAIL.md`](75_C3_REAL_TAIL.md).  Simultaneous coordinate vanishing
remains a separate open obligation.
