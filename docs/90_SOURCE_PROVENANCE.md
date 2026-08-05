# Source provenance

## Primary finite-operator sources

The four programs under `laboratory/` are the preserved research instruments:

1. real float64 all-base scanner;
2. arbitrary-precision Ridder ladder;
3. dyadic asymptotic tail corrector;
4. all-base Green/Bessel radial atlas.

They were imported without changing the finite state, brackets, camera geometry,
or score.  Only filenames and the precision ladder's default sibling-script name
were normalized for repository use.

## Evidence ledgers

`results/c3/` contains direct outputs from the precision and cutoff runs discussed
in the research log.  The large multibase scan transcript is preserved as an XZ-compressed text ledger
under `results/scans/`.

## Contextual notes

The C2 defects and Markov documents are kept under `docs/notes/` because they may
inform future finite-cutoff and return-operator work.  They are not treated as a
derivation of the C3 asymptotic corrector.  Their own evidence-level and firewall
language remains authoritative within those documents.

## Formalization boundary

The Lean root formalizes the finite camera bookkeeping, the real-plane operator,
the faithful quadratic zero detector, and certificate interfaces.  It does not
import numerical JSON values as axioms.
