# Finite native carry operator

## State

For \(n\ge1\),

\[
\psi_t(n)=n^{-1/2}R(-t\log n)e_1.
\]

The coordinate field is \(\mathbb R^2\), and every construction and evaluation
in this repository remains in that real plane.

## Centered bracket

For center \(c\) and radius \(r\),

\[
B_{c,r}(t)=\psi_t(c-r)-2\psi_t(c)+\psi_t(c+r).
\]

## Resultant

With seed coordinates and bracket coordinates indexed by \(e\),

\[
R_{b,M}(t)=\sum_e z_e(t).
\]

The exact finite zero predicate is \(R_{b,M}(t)=0\).  The scanner's score
normalizes \(\lVert R\rVert^2\) by emitted coordinate count and total emitted
energy; the normalization preserves an exact zero whenever its denominator is
positive.

## Objective used by the precision ladder

The fast `resultant` objective minimizes

\[
F_M(t)=\lVert R_{b,M}(t)\rVert^2.
\]

Its stationary numerator is

\[
\frac12F_M'(t)=R_{b,M}(t)\cdot R'_{b,M}(t).
\]

Ridder's bracketed method locates a sign-changing root of this scalar stationary
condition.  Thus the reported number is a finite stationary minimum unless the
resultant coordinates themselves are additionally enclosed at zero.
