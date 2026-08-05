# Oriented C3 tail

## Purpose

The first C3 tail certificate bounds every omitted bracket by its Euclidean
norm.  That proves convergence and a unique limiting minimum, but discards the
rotation responsible for the strong cancellation seen in the real plane.

The oriented certificate keeps both coordinates and encloses only the omitted
sixth-derivative remainder.

## Boundary-jet formula

Let

\[
F_k(x,t)=\partial_t^k\psi_t(x),
\qquad k=0,1,2,
\]

and let \(C=3(M+1)\) be the first omitted C3 center.  Symmetric Taylor
expansion gives

\[
\Delta_1^2F_k(c)
=F_k''(c)+\frac1{12}F_k^{(4)}(c)+\rho_k(c),
\]

with

\[
\lVert\rho_k(c)\rVert
\le\frac1{360}
\sup_{x\in[c-1,c+1]}\lVert F_k^{(6)}(x)\rVert.
\]

Real Euler–Maclaurin bounds, applied separately to the two displayed lattice
sums, produce the oriented boundary jet

\[
\begin{aligned}
T_M^{(k)}(t)={}&
-\frac13F_k'(C,t)
+\frac12F_k''(C,t)
-\frac5{18}F_k'''(C,t)\\
&+\frac1{24}F_k^{(4)}(C,t)
+\frac1{60}F_k^{(5)}(C,t)
+E_{k,M}(t).
\end{aligned}
\]

No fitted coefficient enters this formula.  Every boundary derivative is
computed from the native real rotation.

## Sixth-derivative remainder

For \(|t|\le T\), define

\[
P_j(T)=\prod_{r=0}^{j-1}
\sqrt{T^2+(r+\tfrac12)^2}.
\]

Leibniz differentiation of the factors \(1\), \(-\log x\), and
\((\log x)^2\) gives explicit majorants for
\(\partial_x^6\partial_t^k\psi_t(x)\).  Their integrals are combinations of

\[
\int_A^\infty x^{-13/2}(\log x)^q\,dx,
\qquad q=0,1,2.
\]

The combined real remainder is

\[
\lVert E_{k,M}\rVert
\le\frac7{120}\int_C^\infty
 \lVert\partial_x^6F_k(x)\rVert\,dx
+\frac1{1080}\int_{3M-1}^\infty
 \lVert\partial_x^6F_k(x)\rVert\,dx.
\]

At \(M=16\,364\), on the certified interval, the resulting remainder caps are

\[
\begin{aligned}
\lVert E_{0,M}\rVert &<1.071832\times10^{-16},\\
\lVert E_{1,M}\rVert &<1.184362\times10^{-15},\\
\lVert E_{2,M}\rVert &<1.309052\times10^{-14}.
\end{aligned}
\]

The leading tail is no longer treated as uncertainty: it is the explicitly
enclosed two-coordinate boundary jet.

## Certified limiting minimum

The ledger
[`certificates/c3/c3_m16364_oriented_limit_minimum.json`](../certificates/c3/c3_m16364_oriented_limit_minimum.json)
uses the interval centered at

\[
92.491899270558484305857220904387963299360620986373
\]

with decimal radius \(5\times10^{-16}\).  Directed rounding proves

1. \(H_\infty<0\) at the left endpoint;
2. \(H_\infty>0\) at the right endpoint;
3. \(H'_\infty>0\) throughout the interval;
4. the first coordinate of \(R'_\infty\) is strictly negative there.

Thus the interval contains exactly one stationary point of the limiting
energy, it is a strict minimum, and its velocity is nonzero.

Using the center enclosure and a mean-value bound up to that stationary point,
the certificate also proves

\[
\lVert R_\infty(t_*)\rVert
<2.510236\times10^{-15},
\]

and consequently

\[
|\det(R_\infty(t_*),R'_\infty(t_*))|
<1.161961\times10^{-14}.
\]

These are rigorous positive upper bounds.  They do not establish equality to
zero.  The certificate therefore retains
`limiting_vector_zero_certified: false`.

## Cutoff reinforcement

The same theorem was independently rebuilt with \(M=65\,536\).  The refined
ledger
[`certificates/c3/c3_m65536_oriented_limit_minimum.json`](../certificates/c3/c3_m65536_oriented_limit_minimum.json)
certifies the limiting minimum in the interval centered at

\[
92.4918992705584842962597155468528884029675705307387422925724637743157336566238081857676702969669843062152109549549719933
\]

with radius \(2\times10^{-20}\).  It proves

\[
\lVert R_\infty(t_*)\rVert
<1.661316\times10^{-19},
\]

and

\[
|\det(R_\infty(t_*),R'_\infty(t_*))|
<7.690048\times10^{-19}.
\]

The sixth-derivative remainder for the resultant contracts from approximately
\(1.07\times10^{-16}\) at \(M=16\,364\) to
\(5.20\times10^{-20}\) at \(M=65\,536\).  Both decompositions enclose the same
limiting operator and certify compatible nested locations.  This is a
cutoff-independent reinforcement, not an exact-zero claim.

A third ledger at \(M=131\,072\),
[`c3_m131072_oriented_limit_minimum.json`](../certificates/c3/c3_m131072_oriented_limit_minimum.json),
certifies the same stationary point in a radius \(4\times10^{-22}\) interval.
It proves

\[
\lVert R_\infty(t_*)\rVert<3.477616\times10^{-21}
\]

and contracts the sharp oriented resultant remainder to
\(1.149051\times10^{-21}\).  The separate residual decomposition in
[`78_C3_UNIFORM_RESIDUAL.md`](78_C3_UNIFORM_RESIDUAL.md) identifies which
part of this contraction already has a uniform proof.

## Reverification

```bash
python3 scripts/verify_c3_oriented_tail_certificate.py \
  certificates/c3/c3_m16364_oriented_limit_minimum.json --recompute

python3 scripts/verify_c3_oriented_tail_certificate.py \
  certificates/c3/c3_m65536_oriented_limit_minimum.json --recompute

python3 scripts/verify_c3_oriented_tail_certificate.py \
  certificates/c3/c3_m131072_oriented_limit_minimum.json --recompute
```

The verifier rebuilds the fifth-order boundary jet, all sixth-derivative
integrals, the limiting stationary enclosures, source hashes, and claims.
