# Numerical methods

## Discovery: float64 CPU/CUDA

The discovery scanner evaluates many heights at once.  Its CUDA kernel keeps
all accumulation in IEEE binary64 registers and does not allocate a full
`grid × cutoff` phase matrix.  CPU candidates are re-evaluated canonically in
NumPy float64.

The GPU is therefore a wide-field locator.  At heights near 92, adjacent
float64 values are separated by approximately \(1.4210854715\times10^{-14}\).

## Microscopy: arbitrary precision

The precision ladder reconstructs the same finite state and brackets with
`mpmath`.  It does not call Newton.  It brackets the stationary numerator and
uses Ridder's method, expanding the bracket only when necessary.

## Separation of errors

- **Arithmetic agreement** compares successive decimal-precision stages at a
  fixed cutoff.
- **Cutoff agreement** compares minima of different finite operators.
- **Model agreement** compares an asymptotic prediction against later holdout
  cutoffs.

The practical interpretation must never take more digits than the weakest
relevant layer supports.
