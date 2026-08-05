#!/usr/bin/env python3
"""Verify the oriented real C3 tail and limiting-minimum ledger."""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Sequence


SCHEMA = "org.native-carry.real-oriented-tail-limit-minimum/v1"


def _integer_bounds(enclosure: dict[str, Any]) -> tuple[int, int, int]:
    midpoint = int(enclosure["midpoint_integer"])
    radius = int(enclosure["radius_integer"])
    exponent = int(enclosure["exponent10"])
    if radius < 0:
        raise AssertionError("an enclosure radius cannot be negative")
    return midpoint - radius, midpoint + radius, exponent


def _strictly_negative(enclosure: dict[str, Any]) -> bool:
    _, upper, _ = _integer_bounds(enclosure)
    return upper < 0


def _strictly_positive(enclosure: dict[str, Any]) -> bool:
    lower, _, _ = _integer_bounds(enclosure)
    return lower > 0


def _verify_upper_record(record: dict[str, Any]) -> None:
    _, enclosure_upper, exponent = _integer_bounds(record["formula_enclosure"])
    upper = record["certified_upper"]
    if int(upper["significand_integer"]) != enclosure_upper:
        raise AssertionError("certified upper is not its enclosure endpoint")
    if int(upper["exponent10"]) != exponent:
        raise AssertionError("certified upper exponent mismatch")
    if enclosure_upper <= 0:
        raise AssertionError("certified upper must be positive")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_oriented_certificate(
    path: Path,
    source_root: Path | None = None,
    *,
    recompute: bool = False,
) -> dict[str, Any]:
    certificate = json.loads(path.read_text(encoding="utf-8"))
    if certificate["schema"] != SCHEMA:
        raise AssertionError("unsupported oriented C3 certificate schema")
    if certificate["scope"]["coordinate_field"] != "R^2":
        raise AssertionError("certificate is not in the native real plane")
    if certificate["scope"]["theory_scope"] != "native_real_operator_only":
        raise AssertionError("certificate left the native real theory")
    if certificate["scope"]["oriented_tail"] is not True:
        raise AssertionError("certificate discarded the oriented tail")
    if int(certificate["operator"]["camera"]) != 3:
        raise AssertionError("the oriented theorem is specific to C3")

    domain = certificate["domain"]
    lower = Decimal(domain["exact_decimal_lower"])
    upper = Decimal(domain["exact_decimal_upper"])
    center = Decimal(domain["requested_center"])
    radius = Decimal(domain["requested_radius"])
    with localcontext() as decimal_context:
        decimal_context.prec = max(
            len(domain["exact_decimal_lower"]), len(domain["exact_decimal_upper"])
        ) + 40
        expected_lower = center - radius
        expected_upper = center + radius
    if not (lower == expected_lower and upper == expected_upper and lower < upper):
        raise AssertionError("oriented certificate domain is inconsistent")
    if Decimal(domain["time_abs_upper"]) != max(
        lower.copy_abs(), upper.copy_abs()
    ):
        raise AssertionError("time_abs_upper is not the exact domain cap")

    theorem = certificate["oriented_tail_theorem"]
    expected_preconditions = {
        "cutoff_at_least_two": int(certificate["operator"]["finite_cutoff"]) >= 2,
        "sixth_derivative_majorants_integrable": True,
        "majorants_decrease_on_tail_domain": True,
    }
    if theorem["verified_preconditions"] != expected_preconditions:
        raise AssertionError("oriented-tail preconditions are inconsistent")
    if not all(expected_preconditions.values()):
        raise AssertionError("oriented-tail preconditions were not established")
    if int(theorem["boundary_center"]) != 3 * (
        int(certificate["operator"]["finite_cutoff"]) + 1
    ):
        raise AssertionError("incorrect first omitted C3 center")

    tail_enclosures = certificate["tail_enclosures_on_domain"]
    expected_orders = {
        "resultant": 0,
        "first_time_derivative": 1,
        "second_time_derivative": 2,
    }
    for key, order in expected_orders.items():
        tail = tail_enclosures[key]
        if int(tail["time_derivative_order"]) != order:
            raise AssertionError("oriented tail derivative order mismatch")
        _verify_upper_record(tail["remainder_norm_upper"])

    enclosures = certificate["limit_enclosures"]
    computed_conditions = {
        "limiting_H_left_strictly_negative": _strictly_negative(
            enclosures["H_at_lower"]
        ),
        "limiting_H_right_strictly_positive": _strictly_positive(
            enclosures["H_at_upper"]
        ),
        "limiting_H_prime_strictly_positive_on_domain": _strictly_positive(
            enclosures["H_prime_on_domain"]
        ),
        "limiting_first_derivative_x_strictly_negative_on_domain": (
            _strictly_negative(enclosures["first_derivative_x_on_domain"])
        ),
    }
    if computed_conditions != certificate["verified_conditions"]:
        raise AssertionError("stored oriented conditions do not match enclosures")

    unique_minimum = all(
        computed_conditions[key]
        for key in (
            "limiting_H_left_strictly_negative",
            "limiting_H_right_strictly_positive",
            "limiting_H_prime_strictly_positive_on_domain",
        )
    )
    derivative_nonzero = computed_conditions[
        "limiting_first_derivative_x_strictly_negative_on_domain"
    ]
    expected_claims = {
        "oriented_C3_tail_enclosed_through_second_time_derivative": True,
        "unique_limiting_stationary_point_in_domain": unique_minimum,
        "strict_limiting_minimum_in_domain": unique_minimum,
        "limiting_first_derivative_nonzero_on_domain": derivative_nonzero,
        "stationary_resultant_has_certified_small_norm": True,
        "stationary_determinant_has_certified_small_absolute_value": True,
        "limiting_vector_zero_certified": False,
    }
    if certificate["claims"] != expected_claims:
        raise AssertionError("oriented claims do not follow from the ledger")
    if not (unique_minimum and derivative_nonzero):
        raise AssertionError("oriented ledger does not establish its minimum")

    bounds = certificate["stationary_point_bounds"]
    _verify_upper_record(bounds["resultant_norm_upper"])
    _verify_upper_record(bounds["determinant_abs_upper"])

    root = source_root or path.resolve().parents[2]
    for source in certificate["source_hashes"].values():
        source_path = root / source["path"]
        if _sha256(source_path) != source["sha256"]:
            raise AssertionError(f"source digest mismatch: {source_path}")

    if recompute:
        import sys

        root_text = str(root)
        if root_text not in sys.path:
            sys.path.insert(0, root_text)
        from certification.c3_oriented_tail import (
            certify_oriented_c3_limit_minimum,
        )

        rebuilt = certify_oriented_c3_limit_minimum(
            cutoff=int(certificate["operator"]["finite_cutoff"]),
            center=domain["requested_center"],
            radius=domain["requested_radius"],
            dps=int(certificate["arithmetic"]["decimal_digits"]),
            source_root=root,
        )
        if rebuilt != certificate:
            differing = sorted(
                key for key in rebuilt if rebuilt.get(key) != certificate.get(key)
            )
            raise AssertionError(f"recomputed oriented mismatch: {differing}")

    return certificate


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument(
        "--recompute",
        action="store_true",
        help="repeat every directed-rounding oriented-tail evaluation",
    )
    args = parser.parse_args(argv)
    certificate = verify_oriented_certificate(
        args.certificate, recompute=args.recompute
    )
    print(f"verified={args.certificate}")
    print(f"recomputed={str(args.recompute).lower()}")
    for key, value in certificate["claims"].items():
        print(f"{key}={str(value).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
