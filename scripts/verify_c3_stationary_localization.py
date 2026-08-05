#!/usr/bin/env python3
"""Verify the C3 stationary-localization radius ledger."""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any, Sequence


SCHEMA = "org.native-carry.c3-stationary-localization/v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _integer_bounds(enclosure: dict[str, Any]) -> tuple[int, int, int]:
    midpoint = int(enclosure["midpoint_integer"])
    radius = int(enclosure["radius_integer"])
    exponent = int(enclosure["exponent10"])
    if radius < 0:
        raise AssertionError("an enclosure radius cannot be negative")
    return midpoint - radius, midpoint + radius, exponent


def _verify_upper_record(record: dict[str, Any]) -> None:
    _, enclosure_upper, exponent = _integer_bounds(record["formula_enclosure"])
    upper = record["certified_upper"]
    if int(upper["significand_integer"]) != enclosure_upper:
        raise AssertionError("certified upper is not its enclosure endpoint")
    if int(upper["exponent10"]) != exponent or enclosure_upper <= 0:
        raise AssertionError("invalid positive upper record")


def _verify_lower_record(record: dict[str, Any]) -> None:
    enclosure_lower, _, exponent = _integer_bounds(record["formula_enclosure"])
    lower = record["certified_lower"]
    if int(lower["significand_integer"]) != enclosure_lower:
        raise AssertionError("certified lower is not its enclosure endpoint")
    if int(lower["exponent10"]) != exponent or enclosure_lower <= 0:
        raise AssertionError("invalid positive lower record")


def _decimal_upper(record: dict[str, Any]) -> Decimal:
    upper = record["certified_upper"]
    return Decimal(upper["significand_integer"]) * (
        Decimal(10) ** int(upper["exponent10"])
    )


def _decimal_lower(record: dict[str, Any]) -> Decimal:
    lower = record["certified_lower"]
    return Decimal(lower["significand_integer"]) * (
        Decimal(10) ** int(lower["exponent10"])
    )


def verify_stationary_localization(
    path: Path,
    source_root: Path | None = None,
    *,
    recompute: bool = False,
) -> dict[str, Any]:
    getcontext().prec = 1000
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload["schema"] != SCHEMA:
        raise AssertionError("unsupported stationary-localization schema")
    scope = payload["scope"]
    if scope["coordinate_field"] != "R^2":
        raise AssertionError("stationary localization is not real-plane native")
    if scope["theory_scope"] != "native_real_operator_only":
        raise AssertionError("stationary localization left the native theory")
    if int(scope["camera"]) != 3:
        raise AssertionError("stationary localization is not for C3")
    if scope["shared_limiting_stationary_point"] is not True:
        raise AssertionError("the ledger does not fix one limiting stationary point")

    anchor = payload["anchor"]
    slope_record = anchor["limiting_stationary_slope_lower_m"]
    _verify_lower_record(slope_record)
    slope = _decimal_lower(slope_record)
    for key in (
        "limiting_resultant_norm_upper",
        "limiting_velocity_norm_upper",
        "corrected_resultant_uniform_cap",
        "corrected_velocity_uniform_cap",
    ):
        _verify_upper_record(anchor[key])
    corrected_resultant_cap = _decimal_upper(
        anchor["corrected_resultant_uniform_cap"]
    )
    corrected_velocity_cap = _decimal_upper(
        anchor["corrected_velocity_uniform_cap"]
    )

    perturbation = payload["stationary_perturbation"]
    _verify_upper_record(
        perturbation["resultant_remainder_polynomial_constant"]
    )
    _verify_upper_record(
        perturbation["velocity_remainder_polynomial_constant"]
    )
    if perturbation["corrected_equation"] != "h_M = A_M dot B_M":
        raise AssertionError("unexpected corrected stationary equation")

    record_keys = (
        "corrected_resultant_norm_Q_M_upper",
        "corrected_velocity_norm_upper",
        "corrected_stationary_center_residual_upper",
        "sharp_resultant_remainder_eta0_upper",
        "sharp_velocity_remainder_eta1_upper",
        "stationary_tail_perturbation_upper",
        "total_stationary_error_at_center_upper",
        "derived_localization_radius_upper",
        "certified_radius_upper",
        "resultant_remainder_polynomial_witness_upper",
        "velocity_remainder_polynomial_witness_upper",
        "ideal_root_stationary_error_polynomial_witness_upper",
        "ideal_root_localization_radius_polynomial_witness_upper",
    )
    entries = payload["entries"]
    if len(entries) < 2:
        raise AssertionError("stationary localization has too few entries")
    cutoffs = [int(entry["cutoff"]) for entry in entries]
    if cutoffs != sorted(set(cutoffs)) or cutoffs[0] < 2:
        raise AssertionError("stationary-localization cutoffs are invalid")

    root = source_root or path.resolve().parents[2]
    derived_radii: list[Decimal] = []
    ideal_radii: list[Decimal] = []
    for entry in entries:
        for key in record_keys:
            _verify_upper_record(entry[key])
        certificate_path = root / entry["oriented_certificate_path"]
        if _sha256(certificate_path) != entry["oriented_certificate_sha256"]:
            raise AssertionError(f"oriented certificate digest mismatch: {certificate_path}")

        q_value = _decimal_upper(
            entry["corrected_resultant_norm_Q_M_upper"]
        )
        velocity = _decimal_upper(entry["corrected_velocity_norm_upper"])
        sigma = _decimal_upper(
            entry["corrected_stationary_center_residual_upper"]
        )
        eta0 = _decimal_upper(
            entry["sharp_resultant_remainder_eta0_upper"]
        )
        eta1 = _decimal_upper(entry["sharp_velocity_remainder_eta1_upper"])
        tail = _decimal_upper(entry["stationary_tail_perturbation_upper"])
        total = _decimal_upper(
            entry["total_stationary_error_at_center_upper"]
        )
        radius = _decimal_upper(entry["derived_localization_radius_upper"])
        certified_radius = Decimal(entry["certified_radius"])
        eta0_polynomial = _decimal_upper(
            entry["resultant_remainder_polynomial_witness_upper"]
        )
        eta1_polynomial = _decimal_upper(
            entry["velocity_remainder_polynomial_witness_upper"]
        )
        ideal_error = _decimal_upper(
            entry["ideal_root_stationary_error_polynomial_witness_upper"]
        )
        ideal_radius = _decimal_upper(
            entry["ideal_root_localization_radius_polynomial_witness_upper"]
        )

        if tail < q_value * eta1 + velocity * eta0 + eta0 * eta1:
            raise AssertionError("stationary tail perturbation is too small")
        if total < sigma + tail:
            raise AssertionError("total stationary error is too small")
        if radius * slope < total:
            raise AssertionError("derived localization radius is too small")
        if radius > certified_radius:
            raise AssertionError("derived radius leaves the oriented interval")
        if entry["derived_radius_within_certified_radius"] is not True:
            raise AssertionError("stored radius-containment claim is false")
        if eta0 > eta0_polynomial or eta1 > eta1_polynomial:
            raise AssertionError("a polynomial remainder witness is too small")
        ideal_formula = (
            corrected_resultant_cap * eta1_polynomial
            + corrected_velocity_cap * eta0_polynomial
            + eta0_polynomial * eta1_polynomial
        )
        if ideal_error < ideal_formula:
            raise AssertionError("ideal-root stationary witness is too small")
        if ideal_radius * slope < ideal_error:
            raise AssertionError("ideal-root localization witness is too small")
        derived_radii.append(radius)
        ideal_radii.append(ideal_radius)

    derived_contract = all(
        current < previous
        for previous, current in zip(derived_radii, derived_radii[1:])
    )
    ideal_contract = all(
        current < previous
        for previous, current in zip(ideal_radii, ideal_radii[1:])
    )
    expected_claims = {
        "finite_stationary_localization_radii_certified": True,
        "all_derived_radii_fit_inside_oriented_certificates": True,
        "derived_radii_strictly_contract_on_finite_entries": derived_contract,
        "resultant_and_velocity_remainders_have_polynomial_witnesses": True,
        "second_derivative_remainder_has_polynomial_witness": True,
        "ideal_corrected_root_radius_witness_tends_to_zero": True,
        "corrected_stationary_root_family_constructed_from_threshold": True,
        "limiting_vector_zero_certified": False,
    }
    if payload["claims"] != expected_claims:
        raise AssertionError("stationary-localization claims exceed the ledger")
    if not (derived_contract and ideal_contract):
        raise AssertionError("stationary-localization radii do not contract")

    lean_bridge = payload["lean_bridge"]
    if lean_bridge != {
        "contract": "StationaryLocalizationCertificate",
        "radius_theorem": (
            "StationaryLocalizationCertificate.witness_distance_le_radius"
        ),
        "limit_theorem": "StationaryLocalizationCertificate.radius_tendsToZero",
    }:
        raise AssertionError("unexpected stationary-localization Lean bridge")

    for source in payload["source_hashes"].values():
        source_path = root / source["path"]
        if _sha256(source_path) != source["sha256"]:
            raise AssertionError(f"source digest mismatch: {source_path}")

    if recompute:
        import sys

        root_text = str(root)
        if root_text not in sys.path:
            sys.path.insert(0, root_text)
        from certification.c3_stationary_localization import (
            build_stationary_localization,
        )

        rebuilt = build_stationary_localization(
            [root / entry["oriented_certificate_path"] for entry in entries],
            source_root=root,
        )
        if rebuilt != payload:
            differing = sorted(
                key for key in rebuilt if rebuilt.get(key) != payload.get(key)
            )
            raise AssertionError(
                f"recomputed stationary-localization mismatch: {differing}"
            )

    return payload


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args(argv)
    payload = verify_stationary_localization(
        args.ledger, recompute=args.recompute
    )
    print(f"verified={args.ledger}")
    print(f"recomputed={str(args.recompute).lower()}")
    for key, value in payload["claims"].items():
        print(f"{key}={str(value).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
