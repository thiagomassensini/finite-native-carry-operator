# Finite Native Carry Operator

Finite Native Carry Operator is an audit-oriented laboratory and Lean 4
formalization for the **finite**, mass-built native carry operator in
\(\mathbb R^2\).  It collects the exact camera geometry, the CPU/CUDA discovery
scanner, arbitrary-precision localization without Newton, cutoff-convergence
ledgers, a dyadic asymptotic tail corrector, and directed-rounding real interval
certificates.

- Release target: `v0.1.0`
- Python: 3.12+
- Lean and Mathlib line: `v4.32.0`
- License: Apache-2.0

## Semantic contract

### 1. The native state is fixed before the camera

For every positive integer \(n\),

\[
\psi_t(n)=n^{-1/2}
  \bigl(\cos(-t\log n),\sin(-t\log n)\bigr).
\]

The amplitude \(n^{-1/2}\) is the quadratic root of the inverse-integer carry
mass.  The finite operator does not fit or insert another weight after the
state has been built.

### 2. The finite operator is a real-plane sum

For emitted coordinates \(z_e\in\mathbb R^2\),

\[
R_{b,M}(t)=\sum_e z_e.
\]

The centered bracket is

\[
\Delta_r^2\psi_t(c)=
\psi_t(c-r)-2\psi_t(c)+\psi_t(c+r).
\]

There is no post-bracket calibration matrix.  The raw finite zero equation is
\(R_{b,M}(t)=(0,0)\).

### 3. Camera geometry is explicit

- **C2:** seed \(\psi_t(1)\), centers \(4,8,\ldots,4M\), radius \(r=1\).
- **Natural camera \(b\ge3\):** seeds \(1,\ldots,\lfloor b/2\rfloor\), centers
  \(b,2b,\ldots,Mb\), and every radius
  \(1,\ldots,\lfloor b/2\rfloor\).
- For even \(b\), the last radius is antipodal.  Its two legs remain distinct.
  Consequently C4 contains both \(r=1\) and \(r=2\); C4 is not C2.

### 4. Score preserves exact zeros but is not a digit counter

The scanner publishes

\[
\operatorname{score}(t)=
\frac{\lVert R_{b,M}(t)\rVert^2}
     {N\sum_e\lVert z_e\rVert^2}.
\]

This normalization does not change exact zeros when the denominator is
nonzero.  A tiny score is not, by itself, a statement about decimal accuracy
in \(t\).

### 5. Three different numerical claims remain separated

1. **Finite arithmetic localization:** where the selected finite objective has
   a stationary minimum at fixed \(M\).
2. **Cutoff convergence:** how that minimum moves as \(M\) grows.
3. **Limit certification:** a theorem requiring directed-rounding enclosures,
   uniqueness, and a rigorous tail bound.

The current working tree contains an interval certificate for the fixed-cutoff
C3 minimum at \(M=16\,364\), uniform real tail bounds through the second time
derivative, and a certified unique strict minimum of the limiting C3 energy.
It does **not** claim an interval-certified limiting vector zero.

## First real interval certificate

The certificate
[`certificates/c3/c3_m16364_raw_minimum.json`](certificates/c3/c3_m16364_raw_minimum.json)
uses real directed-rounding balls to prove, on an interval of decimal radius
\(10^{-40}\):

- a unique stationary point of \(\lVert R_{3,16364}(t)\rVert^2\);
- strict positivity of the stationary derivative throughout the interval;
- a strict finite minimum;
- absence of a finite vector zero throughout the same interval.

The minimum claim and the vector-zero claim are deliberately separate.  The
certificate contains no external comparison and uses no non-real coordinate.

## First real C3 tail certificate

The ledger
[`certificates/c3/c3_m16364_tail_limit_minimum.json`](certificates/c3/c3_m16364_tail_limit_minimum.json)
uses the positive-kernel identity for a centered second difference to prove
uniform bounds for the omitted C3 resultant and its first two time derivatives.
At \(M=16\,364\), on an interval of radius \(5\times10^{-5}\), it proves

- \(\lVert T_M\rVert < 1.748059\times10^{-4}\);
- \(\lVert T'_M\rVert < 2.008491\times10^{-3}\);
- \(\lVert T''_M\rVert < 2.315488\times10^{-2}\).

Combining these bounds with a ten-cell directed-rounding cover proves one and
only one stationary point of the limiting C3 energy in the interval, and that
point is a strict minimum.  Simultaneous vanishing of both limiting resultant
coordinates remains an open and separate certification obligation.  The same
cover proves that the limiting velocity is nonzero, so the real-plane Lagrange
identity reduces the remaining obligation exactly to
\(\det(R_\infty,R'_\infty)=0\) at the unique stationary point.
The norm-only tail bound does not decide that determinant.  This motivates the
oriented refinement below, which preserves the cancellation of the C3 tail.

## Oriented C3 tail certificate

The ledger
[`certificates/c3/c3_m16364_oriented_limit_minimum.json`](certificates/c3/c3_m16364_oriented_limit_minimum.json)
implements that refinement.  A fifth-order real boundary jet retains the two
tail coordinates, while explicit sixth-derivative integrals bound the omitted
remainder.  The respective remainder caps for \(R\), \(R'\), and \(R''\) are
approximately \(1.07\times10^{-16}\), \(1.18\times10^{-15}\), and
\(1.31\times10^{-14}\).

It certifies the unique limiting minimum in an interval of radius
\(5\times10^{-16}\) and proves

\[
\lVert R_\infty(t_*)\rVert<2.510236\times10^{-15}.
\]

This is a rigorous upper bound, not an equality proof.  Exact limiting vector
vanishing remains open.

A cutoff-reinforcement ledger at \(M=65\,536\) contracts the stationary
interval radius to \(2\times10^{-20}\) and the certified resultant bound to
\(1.661316\times10^{-19}\):
[`c3_m65536_oriented_limit_minimum.json`](certificates/c3/c3_m65536_oriented_limit_minimum.json).
A third ledger at \(M=131\,072\) contracts them further to radius
\(4\times10^{-22}\) and resultant bound \(3.477616\times10^{-21}\):
[`c3_m131072_oriented_limit_minimum.json`](certificates/c3/c3_m131072_oriented_limit_minimum.json).

## Finite contraction ladder and exact bridge

The nested ledger
[`c3_oriented_contraction_ladder.json`](certificates/c3/c3_oriented_contraction_ladder.json)
verifies that the three oriented certificates enclose the same limiting
stationary point.  The resultant-norm upper bounds and their exact squared
energy bounds contract strictly at both transitions.

The decomposed ledger
[`c3_uniform_residual_decomposition.json`](certificates/c3/c3_uniform_residual_decomposition.json)
now proves that the analytic tail component tends to zero uniformly for fixed
\(T\).  It isolates the remaining bound as

\[
\lVert R_\infty(t_*)\rVert\le Q_M+\eta_M+Vr_M,
\]

where \(Q_M\) is the corrected-center residual, \(\eta_M\to0\) is the
proved oriented-tail remainder, and \(Vr_M\) is localization inside the
nested interval.  The three computed \(Q_M\) bounds contract strictly, but a
cutoff-uniform proof that \(Q_M\to0\) remains open.

The stationary-localization ledger
[`c3_stationary_localization.json`](certificates/c3/c3_stationary_localization.json)
uses the common slope margin \(H_\infty'>21.426618\) to derive the radius

\[
|t_*-c_M|\le
\frac{|h_M(c_M)|+\text{stationary tail perturbation}}{21.426618}.
\]

All three derived radii fit inside the oriented intervals.  For ideal corrected
roots \(h_M(c_M)=0\), the radius has a proved polynomial limit zero.  A
cutoff-uniform construction of those exact corrected roots is still required
before the localization family can be marked complete.

Lean now contains the exact final bridge: if a cutoff-uniform family bounds
that fixed energy by nonnegative errors tending to zero, then the limiting
real-plane resultant is exactly zero.  The current ledger has three certified
members, not an infinite vanishing family, so it deliberately retains
`limiting_vector_zero_certified: false`.

## Main empirical record

For camera 3, the fixed-cutoff arithmetic ladder at \(M=16\,364\) agrees through
118 significant digits across the requested precision stages.  The cutoff
ladder was then extended through \(M=4\,194\,304\).

The two-layer dyadic tail model reports the empirical corrected center

\[
T_{\mathrm{corr}}=
92.491899270558484305857220904387963299360620986373\ldots
\]

with recommended empirical rounding

\[
92.491899270558484.
\]

Its reported radius \(1.7135\times10^{-16}\) is a **model-consensus radius**, not
an interval proof.  Rolling holdouts and all source ledgers are preserved under
[`results/c3/`](results/c3/).

## Repository map

| Path | Purpose |
|---|---|
| [`laboratory/`](laboratory/) | Finite scanner, precision ladder, asymptotic corrector, and multibase Bessel atlas |
| [`certification/`](certification/) | Real directed-rounding evaluator and certificate generator |
| [`certificates/`](certificates/) | Immutable machine-checkable interval ledgers |
| [`results/c3/`](results/c3/) | Immutable C3 precision, cutoff, holdout, and correction ledgers |
| [`FiniteNativeCarryOperator/`](FiniteNativeCarryOperator/) | Lean definitions and certification contracts |
| [`docs/`](docs/) | Mathematical specification, evidence rules, methods, and roadmap |
| [`audit/`](audit/) | SHA-256 source ledger and result index |
| [`tests/`](tests/) | Reproducibility and schema tests |

## Reproduce the core checks

```bash
python3 -m pip install -r requirements.txt
bash scripts/run_python_audits.sh
lake build --wfail FiniteNativeCarryOperator
```

### Rebuild the real interval certificate

```bash
python3 -m certification.real_interval \
  --camera 3 --cutoff 16364 \
  --center 92.4918997313646729318371647865593003364402220801509593905439244254417537329198192691432369926151160829334558221387136109668994565316748777425877647605599898366177116598 \
  --radius 1e-40 --dps 120 \
  --output certificates/c3/c3_m16364_raw_minimum.json

python3 scripts/verify_interval_certificate.py \
  certificates/c3/c3_m16364_raw_minimum.json --recompute
```

### Rebuild the C3 tail and limiting-minimum certificate

```bash
python3 -m certification.c3_tail \
  --cutoff 16364 \
  --center 92.4918997313646729318371647865593003364402220801509593905439244254417537329198192691432369926151160829334558221387136109668994565316748777425877647605599898366177116598 \
  --radius 0.00005 --subdivisions 10 --dps 100 \
  --output certificates/c3/c3_m16364_tail_limit_minimum.json

python3 scripts/verify_c3_tail_certificate.py \
  certificates/c3/c3_m16364_tail_limit_minimum.json --recompute
```

### Rebuild the oriented C3 tail certificate

```bash
python3 -m certification.c3_oriented_tail \
  --cutoff 16364 \
  --center 92.491899270558484305857220904387963299360620986373 \
  --radius 5e-16 --dps 110 \
  --output certificates/c3/c3_m16364_oriented_limit_minimum.json

python3 scripts/verify_c3_oriented_tail_certificate.py \
  certificates/c3/c3_m16364_oriented_limit_minimum.json --recompute
```

### Rebuild the finite C3 contraction ladder

```bash
python3 -m certification.c3_contraction \
  certificates/c3/c3_m16364_oriented_limit_minimum.json \
  certificates/c3/c3_m65536_oriented_limit_minimum.json \
  certificates/c3/c3_m131072_oriented_limit_minimum.json \
  --output certificates/c3/c3_oriented_contraction_ladder.json

python3 scripts/verify_c3_contraction_ladder.py \
  certificates/c3/c3_oriented_contraction_ladder.json --recompute
```

### Rebuild the cutoff-uniform residual decomposition

```bash
python3 -m certification.c3_uniform_residual \
  certificates/c3/c3_m16364_oriented_limit_minimum.json \
  certificates/c3/c3_m65536_oriented_limit_minimum.json \
  certificates/c3/c3_m131072_oriented_limit_minimum.json \
  --output certificates/c3/c3_uniform_residual_decomposition.json

python3 scripts/verify_c3_uniform_residual.py \
  certificates/c3/c3_uniform_residual_decomposition.json --recompute
```

### Rebuild the stationary-localization ledger

```bash
python3 -m certification.c3_stationary_localization \
  certificates/c3/c3_m16364_oriented_limit_minimum.json \
  certificates/c3/c3_m65536_oriented_limit_minimum.json \
  certificates/c3/c3_m131072_oriented_limit_minimum.json \
  --output certificates/c3/c3_stationary_localization.json

python3 scripts/verify_c3_stationary_localization.py \
  certificates/c3/c3_stationary_localization.json --recompute
```

### Float64 CPU/CUDA discovery

```bash
python3 laboratory/native_carry_primitive_real_operator_all_bases.py \
  --cameras 2,3,4,5,6,7,8,9 \
  --tmin 1 --tmax 100 --grid 0.0001 \
  --cutoff 16364 --backend auto \
  --threshold -1 --prominence 0.1 \
  --show score,energy,geometry,radii
```

### Multiprecision localization without Newton

```bash
python3 laboratory/native_carry_precision_ladder.py \
  --camera 3 --cutoff 16364 \
  --t 92.4919 --half-width 0.0001 \
  --dps 40,60,100,160 \
  --objective resultant \
  --cutoff-ladder 2048,4096,8192,16364,32768,65536 \
  --cutoff-dps 60 \
  --json-out c3_precision_max.json
```

### Rebuild the asymptotic audit

```bash
python3 laboratory/native_carry_asymptotic_corrector.py \
  results/c3/native_carry_cutoff_demo_c3.json \
  results/c3/c3_precision_max.json \
  results/c3/c3_cutoff_push_262144.json \
  results/c3/c3_cutoff_holdout_1048576.json \
  results/c3/c3_cutoff_holdout_4194304.json \
  --camera 3 --decay 1.5 \
  --max-linear-layers 3 --canonical-layers 2 \
  --consensus-windows 3 --recent-holdouts 2 \
  --predict-doublings 4 --dps 100 \
  --json-out /tmp/c3_corrector.json \
  --markdown-out /tmp/C3_CORRECTOR.md
```

## Documentation

| Document | Subject |
|---|---|
| [Scope and evidence](docs/00_SCOPE_AND_EVIDENCE.md) | Definitions, exact statements, numerical evidence, and nonclaims |
| [Finite operator](docs/10_FINITE_OPERATOR.md) | State, bracket, resultant, score, and stationary objective |
| [Camera geometry](docs/20_CAMERA_GEOMETRY.md) | C2, natural odd cameras, even antipodal cameras |
| [Numerical methods](docs/30_NUMERICAL_METHODS.md) | CUDA discovery, mpmath, Ridder, and precision separation |
| [Cutoff convergence](docs/40_PRECISION_AND_CUTOFF.md) | Recorded ladders through \(2^{22}\) |
| [Asymptotic corrector](docs/50_ASYMPTOTIC_CORRECTOR.md) | Dyadic rotation/contraction and rolling holdouts |
| [Green/boundary/return](docs/60_GREEN_BOUNDARY_RETURN.md) | Exact finite reconstruction audit and logical scope |
| [Real interval certificate](docs/70_REAL_INTERVAL_CERTIFICATION.md) | Certified C3 finite minimum and finite-zero exclusion |
| [Real C3 tail](docs/75_C3_REAL_TAIL.md) | Uniform derivative tails and certified limiting minimum |
| [Oriented C3 tail](docs/76_C3_ORIENTED_TAIL.md) | Boundary-jet tail enclosure and narrow limiting minimum |
| [C3 contraction criterion](docs/77_C3_CONTRACTION_CRITERION.md) | Nested certificates and the exact cutoff-uniform zero bridge |
| [C3 uniform residual](docs/78_C3_UNIFORM_RESIDUAL.md) | Three-component bound and the proved tail-decay term |
| [C3 stationary localization](docs/79_C3_STATIONARY_LOCALIZATION.md) | Corrected stationary equation, slope margin, and vanishing-radius witness |
| [Interval roadmap](docs/roadmap/INTERVAL_CERTIFICATION.md) | Directed rounding, derivative enclosures, uniqueness, and tail bounds |
| [Reproducibility](docs/80_REPRODUCIBILITY.md) | Commands, environment, ledgers, and release gate |
| [Provenance](docs/90_SOURCE_PROVENANCE.md) | Origin and role of every imported artifact |

## Citation

Versioned citation metadata is in [`CITATION.cff`](CITATION.cff) and
[`.zenodo.json`](.zenodo.json).  A GitHub release is the publication unit; a
Zenodo DOI exists only after Zenodo has ingested that release.
