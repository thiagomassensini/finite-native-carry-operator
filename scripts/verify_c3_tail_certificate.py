#!/usr/bin/env python3
"""Verify the real C3 tail and limiting-minimum certificate ledger."""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Sequence


SCHEMA = "org.native-carry.real-tail-limit-minimum/v1"


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


def _contains_zero(enclosure: dict[str, Any]) -> bool:
    lower, upper, _ = _integer_bounds(enclosure)
    return lower <= 0 <= upper


def _lower_decimal(enclosure: dict[str, Any]) -> Decimal:
    lower, _, exponent = _integer_bounds(enclosure)
    return Decimal(lower) * (Decimal(10) ** exponent)


def _verify_upper_record(record: dict[str, Any]) -> None:
    _, enclosure_upper, exponent = _integer_bounds(record["formula_enclosure"])
    certified = record["certified_upper"]
    if int(certified["significand_integer"]) != enclosure_upper:
        raise AssertionError("certified upper is not the enclosure upper endpoint")
    if int(certified["exponent10"]) != exponent:
        raise AssertionError("certified upper exponent does not match its enclosure")
    if enclosure_upper <= 0:
        raise AssertionError("a certified upper bound must be positive")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_c3_tail_certificate(
    path: Path,
    source_root: Path | None = None,
    *,
    recompute: bool = False,
) -> dict[str, Any]:
    certificate = json.loads(path.read_text(encoding="utf-8"))
    if certificate["schema"] != SCHEMA:
        raise AssertionError("unsupported C3 tail certificate schema")
    if certificate["scope"]["coordinate_field"] != "R^2":
        raise AssertionError("the certificate is not real-plane native")
    if certificate["scope"]["theory_scope"] != "native_real_operator_only":
        raise AssertionError("the certificate left the native real theory")
    if int(certificate["operator"]["camera"]) != 3:
        raise AssertionError("this tail theorem is specific to camera C3")

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
        raise AssertionError("the exact decimal domain is inconsistent")
    expected_time_cap = max(lower.copy_abs(), upper.copy_abs())
    if Decimal(domain["time_abs_upper"]) != expected_time_cap:
        raise AssertionError("time_abs_upper is not the exact domain cap")

    derivation = certificate["tail_derivation"]
    expected_preconditions = {
        "cutoff_at_least_two": int(certificate["operator"]["finite_cutoff"]) >= 2,
        "lower_integration_point_at_least_five": int(
            derivation["lower_integration_point"]
        )
        >= 5,
        "time_domain_is_bounded": True,
        "log_lower_point_strictly_above_four_fifths": (
            _lower_decimal(derivation["log_lower_integration_point"])
            > Decimal(4) / Decimal(5)
        ),
        "majorants_decrease_on_integration_domain": (
            _lower_decimal(derivation["log_lower_integration_point"])
            > Decimal(4) / Decimal(5)
        ),
    }
    if derivation["verified_preconditions"] != expected_preconditions:
        raise AssertionError("stored tail preconditions are inconsistent")
    if not all(expected_preconditions.values()):
        raise AssertionError("tail-majorant preconditions were not established")
    if int(derivation["lower_integration_point"]) != (
        3 * int(certificate["operator"]["finite_cutoff"]) - 1
    ):
        raise AssertionError("incorrect lower integration point")

    for record in certificate["tail_bounds"].values():
        _verify_upper_record(record)

    bridge = certificate["limit_bridge"]
    for side in ("left_endpoint", "right_endpoint"):
        endpoint = bridge[side]
        _verify_upper_record(endpoint["finite_resultant_norm_upper"])
        _verify_upper_record(endpoint["finite_first_derivative_norm_upper"])
        _verify_upper_record(endpoint["H_tail_perturbation_upper"])

    cover = bridge["positive_slope_cover"]
    if len(cover) != int(domain["subdivision_count"]):
        raise AssertionError("subdivision count does not match the positive-slope cover")
    previous_upper = domain["exact_decimal_lower"]
    for index, cell in enumerate(cover):
        if int(cell["index"]) != index:
            raise AssertionError("positive-slope cells are not consecutively indexed")
        if Decimal(cell["exact_decimal_lower"]) != Decimal(previous_upper):
            raise AssertionError("positive-slope cover has a gap or overlap")
        if Decimal(cell["exact_decimal_lower"]) >= Decimal(
            cell["exact_decimal_upper"]
        ):
            raise AssertionError("positive-slope cell is not ordered")
        previous_upper = cell["exact_decimal_upper"]
        for key in (
            "finite_resultant_norm_upper",
            "finite_first_derivative_norm_upper",
            "finite_second_derivative_norm_upper",
            "H_prime_tail_perturbation_upper",
        ):
            _verify_upper_record(cell[key])
    if Decimal(previous_upper) != upper:
        raise AssertionError("positive-slope cover does not reach the domain endpoint")

    computed_conditions = {
        "limiting_H_left_strictly_negative": _strictly_negative(
            bridge["left_endpoint"]["limiting_H_sign_margin"]
        ),
        "limiting_H_right_strictly_positive": _strictly_positive(
            bridge["right_endpoint"]["limiting_H_sign_margin"]
        ),
        "limiting_H_prime_strictly_positive_on_cover": all(
            _strictly_positive(cell["limiting_H_prime_lower_margin"])
            for cell in cover
        ),
        "limiting_first_derivative_x_strictly_negative_on_cover": all(
            _strictly_negative(
                cell["limiting_first_derivative_x_negative_margin"]
            )
            for cell in cover
        ),
    }
    if computed_conditions != certificate["verified_conditions"]:
        raise AssertionError("stored limit conditions do not match exact enclosures")

    limit_minimum = all(
        computed_conditions[key]
        for key in (
            "limiting_H_left_strictly_negative",
            "limiting_H_right_strictly_positive",
            "limiting_H_prime_strictly_positive_on_cover",
        )
    )
    derivative_nonzero = computed_conditions[
        "limiting_first_derivative_x_strictly_negative_on_cover"
    ]
    expected_claims = {
        "uniform_C3_tail_bound_for_resultant": True,
        "uniform_C3_tail_bound_for_first_time_derivative": True,
        "uniform_C3_tail_bound_for_second_time_derivative": True,
        "C3_limit_is_twice_continuously_differentiable_on_domain": True,
        "unique_limiting_stationary_point_in_domain": limit_minimum,
        "strict_limiting_minimum_in_domain": limit_minimum,
        "limiting_first_derivative_nonzero_on_domain": derivative_nonzero,
        "vector_zero_reduced_to_determinant_at_stationary_point": (
            limit_minimum and derivative_nonzero
        ),
        "limiting_vector_zero_certified": False,
    }
    if certificate["claims"] != expected_claims:
        raise AssertionError("stored limit claims do not follow from the ledger")
    if not (limit_minimum and derivative_nonzero):
        raise AssertionError(
            "certificate does not establish the limiting minimum reduction"
        )

    determinant_reference = certificate["vector_zero_reduction"][
        "finite_cutoff_reference_at_requested_center"
    ]
    _verify_upper_record(determinant_reference["K_tail_perturbation_upper"])
    determinant_decided = not _contains_zero(
        determinant_reference["limiting_K_outer_enclosure"]
    )
    if determinant_reference[
        "current_norm_tail_bound_decides_determinant_sign"
    ] != determinant_decided:
        raise AssertionError("stored determinant-decision boundary is inconsistent")
    if determinant_decided:
        raise AssertionError(
            "this ledger unexpectedly decides a determinant sign at its center"
        )

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
        from certification.c3_tail import certify_c3_tail_and_limit_minimum

        rebuilt = certify_c3_tail_and_limit_minimum(
            cutoff=int(certificate["operator"]["finite_cutoff"]),
            center=domain["requested_center"],
            radius=domain["requested_radius"],
            subdivisions=int(domain["subdivision_count"]),
            dps=int(certificate["arithmetic"]["decimal_digits"]),
            source_root=root,
        )
        if rebuilt != certificate:
            differing = sorted(
                key for key in rebuilt if rebuilt.get(key) != certificate.get(key)
            )
            raise AssertionError(f"recomputed certificate mismatch: {differing}")

    return certificate


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument(
        "--recompute",
        action="store_true",
        help="repeat all directed-rounding tail and finite evaluations",
    )
    args = parser.parse_args(argv)
    certificate = verify_c3_tail_certificate(
        args.certificate, recompute=args.recompute
    )
    print(f"verified={args.certificate}")
    print(f"recomputed={str(args.recompute).lower()}")
    for key, value in certificate["claims"].items():
        print(f"{key}={str(value).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
