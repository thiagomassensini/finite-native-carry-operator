# C3 stationary localization from a uniform slope margin

## Corrected stationary equation

Let

\[
A_M=R_M+J_M^{(0)},\qquad
B_M=R_M'+J_M^{(1)},
\]

where the two \(J_M^{(k)}\) are the oriented boundary jets.  The corrected
stationary equation is

\[
h_M=A_M\mathbin{\cdot}B_M.
\]

Write the exact limiting quantities as

\[
R_\infty=A_M+E_{0,M},\qquad
R_\infty'=B_M+E_{1,M},
\]

with \(\lVert E_{0,M}\rVert\le\eta_{0,M}\) and
\(\lVert E_{1,M}\rVert\le\eta_{1,M}\).  Expanding the real scalar product
gives the stationary perturbation bound

\[
|H_\infty-h_M|
\le
\lVert A_M\rVert\eta_{1,M}
+\lVert B_M\rVert\eta_{0,M}
+\eta_{0,M}\eta_{1,M}.
\]

No directional cancellation is assumed in this inequality.

## Localization by the slope margin

The largest certified interval has one directed-rounding lower bound

\[
H_\infty'(t)>m,
\qquad
m=21.4266183401914\ldots
\]

throughout the interval.  If \(t_*\) is its unique stationary point, the
mean-value inequality gives

\[
m|t_*-c_M|
\le |H_\infty(t_*)-H_\infty(c_M)|
=|H_\infty(c_M)|.
\]

Consequently,

\[
\boxed{
|t_*-c_M|
\le
\frac{|h_M(c_M)|+
\lVert A_M(c_M)\rVert\eta_{1,M}+
\lVert B_M(c_M)\rVert\eta_{0,M}+
\eta_{0,M}\eta_{1,M}}{m}.
}
\]

The ledger
[`c3_stationary_localization.json`](../certificates/c3/c3_stationary_localization.json)
evaluates every term separately.

## Finite certified radii

| \(M\) | Certified interval radius | Derived localization radius |
|---:|---:|---:|
| \(16\,364\) | \(5\times10^{-16}\) | \(3.274\times10^{-17}\) |
| \(65\,536\) | \(2\times10^{-20}\) | \(1.124\times10^{-20}\) |
| \(131\,072\) | \(4\times10^{-22}\) | \(2.483\times10^{-22}\) |

Thus every derived radius fits inside its oriented certificate, and the three
derived radii contract strictly.  At \(M=131\,072\), the corrected
stationary-center residual itself is below \(2.252\times10^{-41}\); the
localization radius is dominated by the rigorous tail perturbation.

## Uniform decay for ideal corrected roots

For fixed \(|t|\le T\), both tail components have explicit polynomial
witnesses

\[
\eta_{0,M}\le\frac{C_0(T)}{(M+1)^5},
\qquad
\eta_{1,M}\le\frac{C_1(T)}{(M+1)^5}.
\]

The first bound is the resultant witness from the previous stage.  For the
velocity remainder, its sixth-derivative integral contains at most one
logarithm.  The elementary inequality \(\log x\le\sqrt{x}\) for
\(x\ge1\) converts that integral to the same fifth-power witness.

Choose an ideal corrected center satisfying \(h_M(c_M)=0\).  Uniform caps
\(\lVert A_M\rVert\le\bar A\) and
\(\lVert B_M\rVert\le\bar B\) then give

\[
|t_*-c_M|
\le\frac1m\left[
\frac{\bar A C_1+\bar B C_0}{(M+1)^5}
+\frac{C_0C_1}{(M+1)^{10}}
\right]
\longrightarrow0.
\]

Lean formalizes the division by the positive slope, containment of the
witness, convergence of the radii, and convergence of a fixed velocity cap
times those radii in `StationaryLocalizationCertificate`.

## Remaining boundary

The radius formula and its limit are complete once a cutoff-uniform family of
exact corrected roots \(h_M(c_M)=0\) is instantiated.  The current ledger
certifies three such centers to finite decimal residuals, but does not promote
them to a root family for every cutoff.  It therefore records
`corrected_stationary_root_family_constructed_for_all_cutoffs: false`.

After that family is constructed, localization is discharged and the only
substantive residual obligation is proving
\(Q_M=\lVert A_M(c_M)\rVert\to0\).  Exact limiting vector vanishing is not
claimed at this stage.

## Reverification

```bash
python3 scripts/verify_c3_stationary_localization.py \
  certificates/c3/c3_stationary_localization.json --recompute
```
