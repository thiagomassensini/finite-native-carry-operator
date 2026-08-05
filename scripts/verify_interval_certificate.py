#!/usr/bin/env python3
"""Verify exact-integer claims in a real interval certificate ledger."""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Sequence


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


def _excludes_zero(enclosure: dict[str, Any]) -> bool:
    lower, upper, _ = _integer_bounds(enclosure)
    return upper < 0 or lower > 0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_certificate(
    path: Path,
    source_root: Path | None = None,
    *,
    recompute: bool = False,
) -> dict[str, Any]:
    certificate = json.loads(path.read_text(encoding="utf-8"))
    if certificate["schema"] != "org.native-carry.real-interval-minimum/v1":
        raise AssertionError("unsupported interval certificate schema")
    if certificate["scope"]["coordinate_field"] != "R^2":
        raise AssertionError("the certificate is not a real-plane certificate")
    if certificate["scope"]["theory_scope"] != "finite_native_real_operator_only":
        raise AssertionError("the certificate left the finite native real theory")

    domain = certificate["domain"]
    lower = Decimal(domain["exact_decimal_lower"])
    upper = Decimal(domain["exact_decimal_upper"])
    center = Decimal(domain["requested_center"])
    radius = Decimal(domain["requested_radius"])
    digit_budget = max(
        len(domain["exact_decimal_lower"]), len(domain["exact_decimal_upper"])
    ) + 32
    with localcontext() as decimal_context:
        decimal_context.prec = digit_budget
        expected_lower = center - radius
        expected_upper = center + radius
    if not (lower == expected_lower and upper == expected_upper and lower < upper):
        raise AssertionError("the exact decimal domain is inconsistent")

    enclosures = certificate["enclosures"]
    computed_conditions = {
        "H_lower_strictly_negative": _strictly_negative(enclosures["H_at_lower"]),
        "H_upper_strictly_positive": _strictly_positive(enclosures["H_at_upper"]),
        "H_prime_strictly_positive_on_domain": _strictly_positive(
            enclosures["H_prime_on_domain"]
        ),
        "resultant_x_excludes_zero": _excludes_zero(
            enclosures["resultant_x_on_domain"]
        ),
        "resultant_y_excludes_zero": _excludes_zero(
            enclosures["resultant_y_on_domain"]
        ),
    }
    if computed_conditions != certificate["verified_conditions"]:
        raise AssertionError("stored conditions do not match exact integer enclosures")

    unique_minimum = (
        computed_conditions["H_lower_strictly_negative"]
        and computed_conditions["H_upper_strictly_positive"]
        and computed_conditions["H_prime_strictly_positive_on_domain"]
    )
    zero_absent = (
        computed_conditions["resultant_x_excludes_zero"]
        or computed_conditions["resultant_y_excludes_zero"]
    )
    expected_claims = {
        "unique_stationary_point_in_domain": unique_minimum,
        "strict_finite_minimum_in_domain": unique_minimum,
        "finite_vector_zero_absent_from_domain": zero_absent,
        "limiting_zero_certified": False,
    }
    if expected_claims != certificate["claims"]:
        raise AssertionError("stored claims do not follow from the verified conditions")
    if not unique_minimum:
        raise AssertionError("certificate does not establish a unique finite minimum")

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
        from certification.real_interval import certify_finite_minimum

        rebuilt = certify_finite_minimum(
            camera=int(certificate["operator"]["camera"]),
            cutoff=int(certificate["operator"]["cutoff"]),
            center=domain["requested_center"],
            radius=domain["requested_radius"],
            dps=int(certificate["arithmetic"]["decimal_digits"]),
            source_root=root,
        )
        compared_sections = (
            "arithmetic",
            "claims",
            "domain",
            "enclosures",
            "operator",
            "scope",
            "source_hashes",
            "verified_conditions",
        )
        for section in compared_sections:
            if rebuilt[section] != certificate[section]:
                raise AssertionError(f"recomputed certificate mismatch: {section}")

    return certificate


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument(
        "--recompute",
        action="store_true",
        help="repeat every directed-rounding evaluation before accepting the ledger",
    )
    args = parser.parse_args(argv)
    certificate = verify_certificate(args.certificate, recompute=args.recompute)
    print(f"verified={args.certificate}")
    print(f"recomputed={str(args.recompute).lower()}")
    for key, value in certificate["claims"].items():
        print(f"{key}={str(value).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
