#!/usr/bin/env python3
"""Build a finite contraction ladder from oriented C3 certificates."""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Sequence

from scripts.verify_c3_oriented_tail_certificate import verify_oriented_certificate


SCHEMA = "org.native-carry.real-oriented-contraction-ladder/v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _exact_upper(record: dict[str, Any]) -> Decimal:
    upper = record["certified_upper"]
    return Decimal(upper["significand_integer"]) * (
        Decimal(10) ** int(upper["exponent10"])
    )


def _exact_upper_payload(record: dict[str, Any]) -> dict[str, Any]:
    upper = record["certified_upper"]
    return {
        "significand_integer": upper["significand_integer"],
        "exponent10": int(upper["exponent10"]),
        "exact_value_law": "significand_integer * 10^exponent10",
    }


def _square_upper_payload(upper: dict[str, Any]) -> dict[str, Any]:
    significand = int(upper["significand_integer"])
    exponent = int(upper["exponent10"])
    return {
        "significand_integer": str(significand * significand),
        "exponent10": 2 * exponent,
        "exact_value_law": "significand_integer * 10^exponent10",
    }


def build_contraction_ladder(
    certificate_paths: Sequence[Path],
    *,
    source_root: Path | None = None,
) -> dict[str, Any]:
    if len(certificate_paths) < 2:
        raise ValueError("a contraction ladder needs at least two certificates")
    root = source_root or Path(__file__).resolve().parents[1]

    loaded: list[tuple[Path, dict[str, Any]]] = []
    for path in certificate_paths:
        resolved = path if path.is_absolute() else root / path
        certificate = verify_oriented_certificate(resolved, root, recompute=False)
        loaded.append((resolved, certificate))
    loaded.sort(key=lambda item: int(item[1]["operator"]["finite_cutoff"]))

    limit_objects = {item[1]["scope"]["limit_object"] for item in loaded}
    if len(limit_objects) != 1:
        raise RuntimeError("oriented certificates do not define the same C3 limit")

    entries: list[dict[str, Any]] = []
    for path, certificate in loaded:
        domain = certificate["domain"]
        bounds = certificate["stationary_point_bounds"]
        resultant_remainder = certificate["tail_enclosures_on_domain"]["resultant"][
            "remainder_norm_upper"
        ]
        resultant_upper = _exact_upper_payload(bounds["resultant_norm_upper"])
        entries.append(
            {
                "cutoff": int(certificate["operator"]["finite_cutoff"]),
                "certificate_path": str(path.relative_to(root)),
                "certificate_sha256": _sha256(path),
                "domain": {
                    "exact_decimal_lower": domain["exact_decimal_lower"],
                    "exact_decimal_upper": domain["exact_decimal_upper"],
                    "requested_center": domain["requested_center"],
                    "requested_radius": domain["requested_radius"],
                },
                "resultant_norm_upper_at_shared_stationary_point": resultant_upper,
                "resultant_energy_upper_at_shared_stationary_point": (
                    _square_upper_payload(resultant_upper)
                ),
                "determinant_abs_upper_at_shared_stationary_point": (
                    _exact_upper_payload(bounds["determinant_abs_upper"])
                ),
                "oriented_resultant_tail_remainder_upper": (
                    _exact_upper_payload(resultant_remainder)
                ),
                "certificate_claims": certificate["claims"],
            }
        )

    transitions: list[dict[str, Any]] = []
    for previous, current in zip(entries, entries[1:]):
        previous_lower = Decimal(previous["domain"]["exact_decimal_lower"])
        previous_upper = Decimal(previous["domain"]["exact_decimal_upper"])
        current_lower = Decimal(current["domain"]["exact_decimal_lower"])
        current_upper = Decimal(current["domain"]["exact_decimal_upper"])
        nested = previous_lower <= current_lower <= current_upper <= previous_upper

        previous_resultant = _exact_upper(
            {
                "certified_upper": previous[
                    "resultant_norm_upper_at_shared_stationary_point"
                ]
            }
        )
        current_resultant = _exact_upper(
            {
                "certified_upper": current[
                    "resultant_norm_upper_at_shared_stationary_point"
                ]
            }
        )
        previous_determinant = _exact_upper(
            {
                "certified_upper": previous[
                    "determinant_abs_upper_at_shared_stationary_point"
                ]
            }
        )
        current_determinant = _exact_upper(
            {
                "certified_upper": current[
                    "determinant_abs_upper_at_shared_stationary_point"
                ]
            }
        )
        previous_remainder = _exact_upper(
            {
                "certified_upper": previous[
                    "oriented_resultant_tail_remainder_upper"
                ]
            }
        )
        current_remainder = _exact_upper(
            {
                "certified_upper": current[
                    "oriented_resultant_tail_remainder_upper"
                ]
            }
        )
        if not (
            nested
            and current_resultant < previous_resultant
            and current_determinant < previous_determinant
            and current_remainder < previous_remainder
        ):
            raise RuntimeError("the requested ledgers do not form a strict ladder")

        with localcontext() as decimal_context:
            decimal_context.prec = 80
            resultant_ratio = previous_resultant / current_resultant
            determinant_ratio = previous_determinant / current_determinant
            remainder_ratio = previous_remainder / current_remainder
        transitions.append(
            {
                "from_cutoff": previous["cutoff"],
                "to_cutoff": current["cutoff"],
                "refined_domain_is_subset": nested,
                "same_stationary_point_reason": "the refined interval is contained in the previous interval; both certificates prove existence and uniqueness for the same limiting stationary equation",
                "resultant_bound_strictly_contracts": True,
                "resultant_energy_bound_strictly_contracts": True,
                "determinant_bound_strictly_contracts": True,
                "tail_remainder_strictly_contracts": True,
                "resultant_contraction_factor_display": format(
                    resultant_ratio, ".20E"
                ),
                "determinant_contraction_factor_display": format(
                    determinant_ratio, ".20E"
                ),
                "tail_remainder_contraction_factor_display": format(
                    remainder_ratio, ".20E"
                ),
            }
        )

    source_hashes = {
        "ladder_builder": {
            "path": str(Path(__file__).resolve().relative_to(root)),
            "sha256": _sha256(Path(__file__).resolve()),
        },
        "oriented_certifier": {
            "path": "certification/c3_oriented_tail.py",
            "sha256": _sha256(root / "certification/c3_oriented_tail.py"),
        },
    }
    return {
        "schema": SCHEMA,
        "status": "CERTIFIED_FINITE_ORIENTED_CONTRACTION_LADDER",
        "scope": {
            "coordinate_field": "R^2",
            "theory_scope": "native_real_operator_only",
            "camera": 3,
            "limit_object": next(iter(limit_objects)),
        },
        "entries": entries,
        "transitions": transitions,
        "claims": {
            "finite_nested_contraction_ladder_certified": True,
            "all_entries_bound_the_same_limiting_stationary_point": True,
            "resultant_upper_bounds_strictly_contract_across_entries": True,
            "resultant_energy_upper_bounds_strictly_contract_across_entries": True,
            "determinant_upper_bounds_strictly_contract_across_entries": True,
            "infinite_vanishing_bound_family_certified": False,
            "limiting_vector_zero_certified": False,
        },
        "exact_zero_bridge": {
            "sufficient_condition": "one fixed limiting stationary point with nonnegative resultant-energy upper bounds tending to zero",
            "lean_contract": "VanishingLimitResidualCertificate.witness_zero",
            "current_missing_obligation": "extend the finite ladder to a cutoff-uniform family and prove its resultant-energy upper bound tends to zero",
        },
        "source_hashes": source_hashes,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a finite nested contraction ladder from C3 ledgers"
    )
    parser.add_argument("certificates", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    payload = build_contraction_ladder(args.certificates)
    root = Path(__file__).resolve().parents[1]
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"contraction ladder written to {output}")
    for key, value in payload["claims"].items():
        print(f"{key}={str(value).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
