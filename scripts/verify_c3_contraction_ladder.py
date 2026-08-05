#!/usr/bin/env python3
"""Verify the finite nested C3 oriented-contraction ladder."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence


SCHEMA = "org.native-carry.real-oriented-contraction-ladder/v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_contraction_ladder(
    path: Path,
    source_root: Path | None = None,
    *,
    recompute: bool = False,
) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload["schema"] != SCHEMA:
        raise AssertionError("unsupported contraction-ladder schema")
    if payload["scope"]["coordinate_field"] != "R^2":
        raise AssertionError("contraction ladder is not real-plane native")
    if payload["scope"]["theory_scope"] != "native_real_operator_only":
        raise AssertionError("contraction ladder left the native real theory")
    if int(payload["scope"]["camera"]) != 3:
        raise AssertionError("contraction ladder is not for C3")

    entries = payload["entries"]
    transitions = payload["transitions"]
    if len(entries) < 2 or len(transitions) != len(entries) - 1:
        raise AssertionError("contraction ladder has inconsistent length")
    cutoffs = [int(entry["cutoff"]) for entry in entries]
    if cutoffs != sorted(set(cutoffs)):
        raise AssertionError("contraction cutoffs are not strictly ordered")
    for entry in entries:
        norm_upper = entry[
            "resultant_norm_upper_at_shared_stationary_point"
        ]
        energy_upper = entry[
            "resultant_energy_upper_at_shared_stationary_point"
        ]
        norm_significand = int(norm_upper["significand_integer"])
        norm_exponent = int(norm_upper["exponent10"])
        if int(energy_upper["significand_integer"]) != norm_significand**2:
            raise AssertionError("energy upper is not the exact squared norm upper")
        if int(energy_upper["exponent10"]) != 2 * norm_exponent:
            raise AssertionError("energy-upper decimal exponent is inconsistent")
        if (
            norm_upper["exact_value_law"]
            != "significand_integer * 10^exponent10"
            or energy_upper["exact_value_law"]
            != "significand_integer * 10^exponent10"
        ):
            raise AssertionError("unsupported exact upper-bound encoding")

    expected_claims = {
        "finite_nested_contraction_ladder_certified": True,
        "all_entries_bound_the_same_limiting_stationary_point": True,
        "resultant_upper_bounds_strictly_contract_across_entries": True,
        "resultant_energy_upper_bounds_strictly_contract_across_entries": True,
        "determinant_upper_bounds_strictly_contract_across_entries": True,
        "infinite_vanishing_bound_family_certified": False,
        "limiting_vector_zero_certified": False,
    }
    if payload["claims"] != expected_claims:
        raise AssertionError("contraction claims exceed the finite ledger")
    for transition in transitions:
        for key in (
            "refined_domain_is_subset",
            "resultant_bound_strictly_contracts",
            "resultant_energy_bound_strictly_contracts",
            "determinant_bound_strictly_contracts",
            "tail_remainder_strictly_contracts",
        ):
            if transition[key] is not True:
                raise AssertionError(f"failed contraction relation: {key}")

    root = source_root or path.resolve().parents[2]
    for entry in entries:
        certificate_path = root / entry["certificate_path"]
        if _sha256(certificate_path) != entry["certificate_sha256"]:
            raise AssertionError(
                f"oriented certificate digest mismatch: {certificate_path}"
            )
    for source in payload["source_hashes"].values():
        source_path = root / source["path"]
        if _sha256(source_path) != source["sha256"]:
            raise AssertionError(f"source digest mismatch: {source_path}")

    if recompute:
        import sys

        root_text = str(root)
        if root_text not in sys.path:
            sys.path.insert(0, root_text)
        from certification.c3_contraction import build_contraction_ladder

        rebuilt = build_contraction_ladder(
            [root / entry["certificate_path"] for entry in entries],
            source_root=root,
        )
        if rebuilt != payload:
            differing = sorted(
                key for key in rebuilt if rebuilt.get(key) != payload.get(key)
            )
            raise AssertionError(f"recomputed contraction mismatch: {differing}")

    return payload


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ladder", type=Path)
    parser.add_argument(
        "--recompute",
        action="store_true",
        help="rebuild all finite contraction relations",
    )
    args = parser.parse_args(argv)
    payload = verify_contraction_ladder(args.ladder, recompute=args.recompute)
    print(f"verified={args.ladder}")
    print(f"recomputed={str(args.recompute).lower()}")
    for key, value in payload["claims"].items():
        print(f"{key}={str(value).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
