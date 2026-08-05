# Real C3 tail and limiting minimum

## Exact tail geometry

For camera C3, the seed is \(\psi_t(1)\) and every center has radius one.
After the finite cutoff \(M\), the omitted real-plane vector is therefore

\[
T_M(t)=\sum_{m=M+1}^{\infty}
\bigl(\psi_t(3m-1)-2\psi_t(3m)+\psi_t(3m+1)\bigr).
\]

No fitted corrector enters this identity.

## Positive-kernel reduction

For a twice differentiable real-plane function \(f\),

\[
f(c-1)-2f(c)+f(c+1)
=\int_{-1}^{1}(1-|u|)f''(c+u)\,du.
\]

The tent kernel is nonnegative and has mass one.  Write

\[
\psi_t(x)=x^{-1/2}\operatorname{Rot}(-t\log x)e_1,
\qquad |t|\le T.
\]

Direct differentiation in the real plane gives

\[
\lVert\partial_x\psi_t(x)\rVert
\le A_1x^{-3/2},
\quad
A_1=\sqrt{T^2+\tfrac14},
\]

\[
\lVert\partial_x^2\psi_t(x)\rVert
\le A_2x^{-5/2},
\quad
A_2=\sqrt{(T^2+\tfrac14)(T^2+\tfrac94)}.
\]

Commuting zero, one, or two time derivatives with the bracket and applying the
real triangle inequality produces the decreasing majorants

\[
B_0(x)=A_2x^{-5/2},
\]

\[
B_1(x)=x^{-5/2}\bigl(A_2\log x+2A_1+1\bigr),
\]

\[
B_2(x)=x^{-5/2}\bigl(A_2(\log x)^2+(4A_1+2)\log x+2\bigr).
\]

For \(M\ge2\), these functions decrease throughout the required integration
domain.  The integral comparison is

\[
\sum_{m=M+1}^{\infty}B_k(3m-1)
\le\frac13\int_{3M-1}^{\infty}B_k(x)\,dx.
\]

All three integrals are explicit combinations of

\[
I_q(A)=\int_A^\infty x^{-5/2}(\log x)^q\,dx,
\qquad q=0,1,2.
\]

Their common factor is \(A^{-3/2}\), multiplied by a polynomial of degree at
most two in \(\log A\).  Consequently all three bounds tend to zero with the
cutoff and prove uniform convergence through the second time derivative.

## Certified bounds at M = 16,364

On the interval centered at the recorded finite minimum with radius
\(5\times10^{-5}\), the directed-rounding ledger proves

\[
\lVert T_M\rVert < 1.748059\times10^{-4},
\]

\[
\lVert T'_M\rVert < 2.008491\times10^{-3},
\]

\[
\lVert T''_M\rVert < 2.315488\times10^{-2}.
\]

The full decimal upper endpoints are stored in
[`certificates/c3/c3_m16364_tail_limit_minimum.json`](../certificates/c3/c3_m16364_tail_limit_minimum.json).

## Transport to the limiting energy

Let \(R_\infty=R_M+T_M\) and
\(H=R\mathbin{\cdot}R'\).  The endpoint perturbation obeys

\[
|H_\infty-H_M|
\le \lVert R_M\rVert\varepsilon_1
 +\lVert R'_M\rVert\varepsilon_0
 +\varepsilon_0\varepsilon_1.
\]

For the slope,

\[
\begin{aligned}
|H'_\infty-H'_M|\le{}&
2\lVert R'_M\rVert\varepsilon_1+\varepsilon_1^2
+\lVert R_M\rVert\varepsilon_2\\
&+\varepsilon_0\lVert R''_M\rVert
+\varepsilon_0\varepsilon_2.
\end{aligned}
\]

Arb proves the limiting stationary numerator negative at the left endpoint and
positive at the right endpoint.  A ten-cell interval cover proves
\(H'_\infty>0\) everywhere; the smallest stored lower margin is greater than
\(18.79\).  Hence the limiting C3 energy has exactly one stationary point in
the interval, and that point is a strict minimum.

## Exact scalar gate for the vector zero

Define, entirely in the oriented real plane,

\[
K_\infty(t)=\det(R_\infty(t),R'_\infty(t))
=X_\infty Y'_\infty-Y_\infty X'_\infty.
\]

The two-dimensional Lagrange identity is

\[
H_\infty(t)^2+K_\infty(t)^2
=\lVert R_\infty(t)\rVert^2
 \lVert R'_\infty(t)\rVert^2.
\]

The interval cover proves more than positive stationary slope: the first
coordinate of \(R'_\infty\) is strictly negative on every cell.  Thus
\(R'_\infty\neq0\) throughout the domain.  At the unique stationary point
\(t_*\), where \(H_\infty(t_*)=0\), the identity gives the exact equivalence

\[
R_\infty(t_*)=0
\quad\Longleftrightarrow\quad
K_\infty(t_*)=0.
\]

This reduces the remaining vector question to one native real scalar, but it
does not prove that this determinant vanishes.  The certificate explicitly
records both the successful reduction and
`limiting_vector_zero_certified: false`.

At the recorded \(M=16\,364\) finite minimum, the ledger encloses

\[
K_M=-8.639736582897062\ldots\times10^{-6}.
\]

Applying only the norm tail bounds gives a determinant perturbation cap of
approximately \(8.0952\times10^{-4}\), so the resulting limiting determinant
enclosure contains zero.  This is the exact quantitative obstruction: the
coordinate-blind tail estimate is almost two orders of magnitude too wide even
to decide a sign at that reference point.  The next improvement must retain
the oriented cancellation of the real tail rather than bounding every bracket
only by its norm.

That oriented improvement is now implemented and certified in
[`76_C3_ORIENTED_TAIL.md`](76_C3_ORIENTED_TAIL.md).  It replaces the leading
tail uncertainty by an explicit boundary jet and bounds only a
sixth-derivative remainder.

## Reverification

```bash
python3 scripts/verify_c3_tail_certificate.py \
  certificates/c3/c3_m16364_tail_limit_minimum.json --recompute
```

Recomputation rebuilds every real ball, every perturbation margin, the complete
subinterval cover, and all source digests.
