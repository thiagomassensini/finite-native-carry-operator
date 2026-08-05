# Reproducibility

## Python environment

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
bash scripts/run_python_audits.sh
```

CUDA is optional.  CPU self-tests always run; CUDA parity is checked when CuPy
and a device are available.

## Lean environment

```bash
curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh
lake update
lake build --wfail FiniteNativeCarryOperator
```

The toolchain and Mathlib release line are pinned.  `lake-manifest.json` is
committed after regeneration.

## Immutable numerical ledgers

Each JSON preserves:

- camera and cutoff;
- coordinate count and geometry;
- resultants, norms, and scores;
- arithmetic method and precision;
- whether Newton or interval arithmetic was used;
- script SHA-256 and runtime metadata.

`audit/SOURCE_SHA256.txt` covers every versioned source and evidence file except
itself.  `audit/RESULT_LEDGER_INDEX.json` records the semantic role of the main
ledgers.
