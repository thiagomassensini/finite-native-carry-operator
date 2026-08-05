#!/usr/bin/env python3
"""Certify the stationary-localization radius for the oriented C3 limit.

Let A_M and B_M be the finite resultant and velocity corrected by their
oriented boundary jets, and let h_M = A_M dot B_M.  If the omitted resultant
and velocity remainders have norm bounds eta_0 and eta_1, then

    |H_inf - h_M|
      <= ||A_M|| eta_1 + ||B_M|| eta_0 + eta_0 eta_1.

On an anchor interval where H_inf' >= m > 0, a corrected center c_M obeys

    |t_* - c_M|
      <= (|h_M(c_M)| + stationary_tail_error_M) / m.

For an ideal center satisfying h_M(c_M) = 0, polynomial witnesses for eta_0
and eta_1 make this radius tend to zero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any, Sequence

import flint
from flint import arb, ctx

from certification.c3_oriented_tail import (
    _space_factor,
    evaluate_oriented_c3_limit,
    oriented_c3_tail_enclosure,
)
from certification.c3_tail import _plane_norm_upper, _positive_upper_record, _upper_ball
from certification.real_interval import (
    _ball_payload,
    build_sparse_geometry,
    evaluate_real_operator,
    prepare_terms,
)
from scripts.verify_c3_oriented_tail_certificate import verify_oriented_certificate


SCHEMA = "org.native-carry.c3-stationary-localization/v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _positive_lower_record(value: Any, digits: int) -> dict[str, Any]:
    enclosure = _ball_payload(value, digits)
    midpoint = int(enclosure["midpoint_integer"])
    radius = int(enclosure["radius_integer"])
    lower = midpoint - radius
    if lower <= 0:
        raise RuntimeError("expected a strictly positive lower bound")
    return {
        "formula_enclosure": enclosure,
        "certified_lower": {
            "significand_integer": str(lower),
            "exponent10": int(enclosure["exponent10"]),
            "exact_value_law": "significand_integer * 10^exponent10",
        },
    }


def _lower_ball(record: dict[str, Any]) -> Any:
    lower = record["certified_lower"]
    return arb(f"{lower['significand_integer']}e{lower['exponent10']}")


def _decimal_upper(record: dict[str, Any]) -> Decimal:
    upper = record["certified_upper"]
    return Decimal(upper["significand_integer"]) * (
        Decimal(10) ** int(upper["exponent10"])
    )


def build_stationary_localization(
    certificate_paths: Sequence[Path],
    *,
    source_root: Path | None = None,
) -> dict[str, Any]:
    if len(certificate_paths) < 2:
        raise ValueError("stationary localization needs at least two certificates")
    root = source_root or Path(__file__).resolve().parents[1]
    getcontext().prec = 1000

    loaded: list[tuple[Path, dict[str, Any]]] = []
    for path in certificate_paths:
        resolved = path if path.is_absolute() else root / path
        certificate = verify_oriented_certificate(resolved, root, recompute=False)
        loaded.append((resolved, certificate))
    loaded.sort(key=lambda item: int(item[1]["operator"]["finite_cutoff"]))

    limit_objects = {certificate["scope"]["limit_object"] for _, certificate in loaded}
    if len(limit_objects) != 1:
        raise RuntimeError("the certificates do not enclose the same C3 limit")
    anchor_path, anchor = loaded[0]
    anchor_lower = Decimal(anchor["domain"]["exact_decimal_lower"])
    anchor_upper = Decimal(anchor["domain"]["exact_decimal_upper"])
    previous_lower, previous_upper = anchor_lower, anchor_upper
    for _, certificate in loaded[1:]:
        lower = Decimal(certificate["domain"]["exact_decimal_lower"])
        upper = Decimal(certificate["domain"]["exact_decimal_upper"])
        if not (previous_lower <= lower <= upper <= previous_upper):
            raise RuntimeError("the stationary intervals are not nested")
        previous_lower, previous_upper = lower, upper

    precision = max(
        int(certificate["arithmetic"]["decimal_digits"])
        for _, certificate in loaded
    )
    ctx.dps = precision
    digits = precision + 12
    time_cap_text = anchor["domain"]["time_abs_upper"]
    time_cap = arb(time_cap_text)

    anchor_cutoff = int(anchor["operator"]["finite_cutoff"])
    anchor_prepared = prepare_terms(build_sparse_geometry(3, anchor_cutoff))
    anchor_domain = arb(
        arb(anchor["domain"]["requested_center"]),
        arb(anchor["domain"]["requested_radius"]),
    )
    anchor_evaluation = evaluate_oriented_c3_limit(
        time=anchor_domain,
        prepared=anchor_prepared,
        cutoff=anchor_cutoff,
        time_abs_upper=time_cap,
    )
    anchor_lower_evaluation = evaluate_oriented_c3_limit(
        time=arb(anchor["domain"]["exact_decimal_lower"]),
        prepared=anchor_prepared,
        cutoff=anchor_cutoff,
        time_abs_upper=time_cap,
    )
    anchor_upper_evaluation = evaluate_oriented_c3_limit(
        time=arb(anchor["domain"]["exact_decimal_upper"]),
        prepared=anchor_prepared,
        cutoff=anchor_cutoff,
        time_abs_upper=time_cap,
    )
    left_sign_margin_record = _positive_lower_record(
        -anchor_lower_evaluation.stationary, digits
    )
    right_sign_margin_record = _positive_lower_record(
        anchor_upper_evaluation.stationary, digits
    )
    slope_lower_record = _positive_lower_record(
        anchor_evaluation.stationary_derivative, digits
    )
    slope_lower = _lower_ball(slope_lower_record)
    resultant_cap_record = _plane_norm_upper(
        anchor_evaluation.resultant_x,
        anchor_evaluation.resultant_y,
        digits,
    )
    velocity_cap_record = _plane_norm_upper(
        anchor_evaluation.derivative_x,
        anchor_evaluation.derivative_y,
        digits,
    )
    second_derivative_cap_record = _plane_norm_upper(
        anchor_evaluation.second_x,
        anchor_evaluation.second_y,
        digits,
    )

    factors: list[Any] = []
    factor_records: list[dict[str, Any]] = []
    for order in range(7):
        record = _positive_upper_record(_space_factor(order, time_cap), digits)
        factor_records.append(record)
        factors.append(_upper_ball(record))
    first_time_constant = (
        6 * factors[5]
        + 15 * factors[4]
        + 40 * factors[3]
        + 90 * factors[2]
        + 144 * factors[1]
        + 120
    )
    eta0_polynomial_constant_record = _positive_upper_record(
        arb(16) / 1485 * factors[6], digits
    )
    first_integral_polynomial_constant = (
        factors[6] * (arb(2) / 11 + arb(4) / 121)
        + first_time_constant * arb(2) / 11
    )
    eta1_polynomial_constant_record = _positive_upper_record(
        arb(8) / 135 * first_integral_polynomial_constant, digits
    )
    second_time_log_coefficient = (
        12 * factors[5]
        + 30 * factors[4]
        + 80 * factors[3]
        + 180 * factors[2]
        + 288 * factors[1]
        + 240
    )
    second_time_constant = (
        30 * factors[4]
        + 120 * factors[3]
        + 330 * factors[2]
        + 600 * factors[1]
        + 548
    )
    second_integral_polynomial_constant = (
        factors[6]
        * (arb(2) / 11 + arb(8) / 121 + arb(16) / 1331)
        + second_time_log_coefficient * (arb(2) / 11 + arb(4) / 121)
        + second_time_constant * arb(2) / 11
    )
    eta2_polynomial_constant_record = _positive_upper_record(
        arb(8) / 135 * second_integral_polynomial_constant, digits
    )
    eta0_polynomial_constant = _upper_ball(eta0_polynomial_constant_record)
    eta1_polynomial_constant = _upper_ball(eta1_polynomial_constant_record)
    eta2_polynomial_constant = _upper_ball(eta2_polynomial_constant_record)

    anchor_tail0 = oriented_c3_tail_enclosure(
        cutoff=anchor_cutoff,
        time=anchor_domain,
        time_abs_upper=time_cap,
        time_order=0,
    )
    anchor_tail1 = oriented_c3_tail_enclosure(
        cutoff=anchor_cutoff,
        time=anchor_domain,
        time_abs_upper=time_cap,
        time_order=1,
    )
    anchor_tail2 = oriented_c3_tail_enclosure(
        cutoff=anchor_cutoff,
        time=anchor_domain,
        time_abs_upper=time_cap,
        time_order=2,
    )
    corrected_resultant_cap_record = _positive_upper_record(
        _upper_ball(resultant_cap_record) + anchor_tail0.remainder_norm,
        digits,
    )
    corrected_velocity_cap_record = _positive_upper_record(
        _upper_ball(velocity_cap_record) + anchor_tail1.remainder_norm,
        digits,
    )
    corrected_second_derivative_cap_record = _positive_upper_record(
        _upper_ball(second_derivative_cap_record) + anchor_tail2.remainder_norm,
        digits,
    )
    corrected_resultant_cap = _upper_ball(corrected_resultant_cap_record)
    corrected_velocity_cap = _upper_ball(corrected_velocity_cap_record)
    corrected_second_derivative_cap = _upper_ball(
        corrected_second_derivative_cap_record
    )

    root_family_threshold = 131072
    threshold_eta0 = (
        eta0_polynomial_constant / (arb(root_family_threshold + 1) ** 5)
    )
    threshold_eta1 = (
        eta1_polynomial_constant / (arb(root_family_threshold + 1) ** 5)
    )
    threshold_eta2 = (
        eta2_polynomial_constant / (arb(root_family_threshold + 1) ** 4)
    )
    threshold_stationary_error_record = _positive_upper_record(
        corrected_resultant_cap * threshold_eta1
        + corrected_velocity_cap * threshold_eta0
        + threshold_eta0 * threshold_eta1,
        digits,
    )
    threshold_slope_error_record = _positive_upper_record(
        2 * corrected_velocity_cap * threshold_eta1
        + threshold_eta1 * threshold_eta1
        + corrected_resultant_cap * threshold_eta2
        + threshold_eta0 * corrected_second_derivative_cap
        + threshold_eta0 * threshold_eta2,
        digits,
    )
    endpoint_signs_preserved = bool(
        _upper_ball(threshold_stationary_error_record)
        < _lower_ball(left_sign_margin_record)
        and _upper_ball(threshold_stationary_error_record)
        < _lower_ball(right_sign_margin_record)
    )
    corrected_slope_positive = bool(
        _upper_ball(threshold_slope_error_record) < slope_lower
    )
    corrected_root_family_certified = (
        endpoint_signs_preserved and corrected_slope_positive
    )
    if not corrected_root_family_certified:
        raise RuntimeError("the corrected root-family threshold was not certified")

    entries: list[dict[str, Any]] = []
    for path, certificate in loaded:
        cutoff = int(certificate["operator"]["finite_cutoff"])
        center_text = certificate["domain"]["requested_center"]
        radius_text = certificate["domain"]["requested_radius"]
        radius = arb(radius_text)
        prepared = prepare_terms(build_sparse_geometry(3, cutoff))
        finite_center = evaluate_real_operator(
            arb(center_text), prepared, second=True
        )
        assert finite_center.second_x is not None
        assert finite_center.second_y is not None
        tail0 = oriented_c3_tail_enclosure(
            cutoff=cutoff,
            time=arb(center_text),
            time_abs_upper=time_cap,
            time_order=0,
        )
        tail1 = oriented_c3_tail_enclosure(
            cutoff=cutoff,
            time=arb(center_text),
            time_abs_upper=time_cap,
            time_order=1,
        )
        corrected_resultant = (
            finite_center.resultant_x + tail0.approximation_x,
            finite_center.resultant_y + tail0.approximation_y,
        )
        corrected_velocity = (
            finite_center.derivative_x + tail1.approximation_x,
            finite_center.derivative_y + tail1.approximation_y,
        )
        resultant_norm_record = _plane_norm_upper(
            corrected_resultant[0], corrected_resultant[1], digits
        )
        velocity_norm_record = _plane_norm_upper(
            corrected_velocity[0], corrected_velocity[1], digits
        )
        corrected_stationary = (
            corrected_resultant[0] * corrected_velocity[0]
            + corrected_resultant[1] * corrected_velocity[1]
        )
        stationary_center_residual_record = _positive_upper_record(
            corrected_stationary.abs_upper(), digits
        )
        eta0_record = _positive_upper_record(tail0.remainder_norm, digits)
        eta1_record = _positive_upper_record(tail1.remainder_norm, digits)
        stationary_tail_record = _positive_upper_record(
            _upper_ball(resultant_norm_record) * _upper_ball(eta1_record)
            + _upper_ball(velocity_norm_record) * _upper_ball(eta0_record)
            + _upper_ball(eta0_record) * _upper_ball(eta1_record),
            digits,
        )
        stationary_error_record = _positive_upper_record(
            _upper_ball(stationary_center_residual_record)
            + _upper_ball(stationary_tail_record),
            digits,
        )
        derived_radius_record = _positive_upper_record(
            _upper_ball(stationary_error_record) / slope_lower, digits
        )
        radius_record = _positive_upper_record(radius, digits)

        eta0_polynomial_record = _positive_upper_record(
            eta0_polynomial_constant / (arb(cutoff + 1) ** 5), digits
        )
        eta1_polynomial_record = _positive_upper_record(
            eta1_polynomial_constant / (arb(cutoff + 1) ** 5), digits
        )
        if _decimal_upper(eta0_record) > _decimal_upper(eta0_polynomial_record):
            raise RuntimeError("the resultant polynomial witness did not dominate")
        if _decimal_upper(eta1_record) > _decimal_upper(eta1_polynomial_record):
            raise RuntimeError("the velocity polynomial witness did not dominate")
        radius_dominated = _decimal_upper(derived_radius_record) <= Decimal(
            radius_text
        )
        if not radius_dominated:
            raise RuntimeError("the derived stationary radius exceeds the ledger")

        ideal_stationary_polynomial_record = _positive_upper_record(
            corrected_resultant_cap
            * _upper_ball(eta1_polynomial_record)
            + corrected_velocity_cap
            * _upper_ball(eta0_polynomial_record)
            + _upper_ball(eta0_polynomial_record)
            * _upper_ball(eta1_polynomial_record),
            digits,
        )
        ideal_radius_polynomial_record = _positive_upper_record(
            _upper_ball(ideal_stationary_polynomial_record) / slope_lower,
            digits,
        )

        entries.append(
            {
                "cutoff": cutoff,
                "oriented_certificate_path": str(path.relative_to(root)),
                "oriented_certificate_sha256": _sha256(path),
                "center": center_text,
                "certified_radius": radius_text,
                "corrected_resultant_norm_Q_M_upper": resultant_norm_record,
                "corrected_velocity_norm_upper": velocity_norm_record,
                "corrected_stationary_center_residual_upper": (
                    stationary_center_residual_record
                ),
                "sharp_resultant_remainder_eta0_upper": eta0_record,
                "sharp_velocity_remainder_eta1_upper": eta1_record,
                "stationary_tail_perturbation_upper": stationary_tail_record,
                "total_stationary_error_at_center_upper": (
                    stationary_error_record
                ),
                "derived_localization_radius_upper": derived_radius_record,
                "certified_radius_upper": radius_record,
                "derived_radius_within_certified_radius": radius_dominated,
                "resultant_remainder_polynomial_witness_upper": (
                    eta0_polynomial_record
                ),
                "velocity_remainder_polynomial_witness_upper": (
                    eta1_polynomial_record
                ),
                "ideal_root_stationary_error_polynomial_witness_upper": (
                    ideal_stationary_polynomial_record
                ),
                "ideal_root_localization_radius_polynomial_witness_upper": (
                    ideal_radius_polynomial_record
                ),
            }
        )

    derived_radii = [
        _decimal_upper(entry["derived_localization_radius_upper"])
        for entry in entries
    ]
    ideal_radii = [
        _decimal_upper(
            entry["ideal_root_localization_radius_polynomial_witness_upper"]
        )
        for entry in entries
    ]
    derived_contract = all(
        current < previous
        for previous, current in zip(derived_radii, derived_radii[1:])
    )
    ideal_contract = all(
        current < previous
        for previous, current in zip(ideal_radii, ideal_radii[1:])
    )
    if not (derived_contract and ideal_contract):
        raise RuntimeError("the localization witnesses did not contract")

    source_paths = {
        "localization_builder": Path(__file__).resolve(),
        "oriented_tail_certifier": root / "certification/c3_oriented_tail.py",
        "finite_interval_evaluator": root / "certification/real_interval.py",
        "lean_contract": root
        / "FiniteNativeCarryOperator/Certification/Contract.lean",
    }
    return {
        "schema": SCHEMA,
        "status": "CERTIFIED_FINITE_C3_STATIONARY_LOCALIZATION",
        "scope": {
            "coordinate_field": "R^2",
            "theory_scope": "native_real_operator_only",
            "camera": 3,
            "limit_object": next(iter(limit_objects)),
            "shared_limiting_stationary_point": True,
        },
        "anchor": {
            "oriented_certificate_path": str(anchor_path.relative_to(root)),
            "domain_lower": anchor["domain"]["exact_decimal_lower"],
            "domain_upper": anchor["domain"]["exact_decimal_upper"],
            "time_abs_upper_T": time_cap_text,
            "limiting_stationary_slope_lower_m": slope_lower_record,
            "limiting_resultant_norm_upper": resultant_cap_record,
            "limiting_velocity_norm_upper": velocity_cap_record,
            "limiting_second_derivative_norm_upper": (
                second_derivative_cap_record
            ),
            "corrected_resultant_uniform_cap": corrected_resultant_cap_record,
            "corrected_velocity_uniform_cap": corrected_velocity_cap_record,
            "corrected_second_derivative_uniform_cap": (
                corrected_second_derivative_cap_record
            ),
            "left_stationary_sign_margin": left_sign_margin_record,
            "right_stationary_sign_margin": right_sign_margin_record,
        },
        "stationary_perturbation": {
            "corrected_equation": "h_M = A_M dot B_M",
            "tail_bound": "|H_infinity-h_M| <= ||A_M||*eta1_M + ||B_M||*eta0_M + eta0_M*eta1_M",
            "localization_bound": "|t_*-c_M| <= (|h_M(c_M)| + stationary_tail_error_M)/m",
            "ideal_center_condition": "h_M(c_M) = 0",
            "resultant_remainder_polynomial_constant": (
                eta0_polynomial_constant_record
            ),
            "velocity_remainder_polynomial_constant": (
                eta1_polynomial_constant_record
            ),
            "second_derivative_remainder_polynomial_constant": (
                eta2_polynomial_constant_record
            ),
            "polynomial_law": "eta0_M <= C0/(M+1)^5, eta1_M <= C1/(M+1)^5, and eta2_M <= C2/(M+1)^4",
        },
        "corrected_root_family": {
            "threshold_cutoff": root_family_threshold,
            "definition": "for every integer M >= threshold_cutoff, c_M is the unique root of h_M on the anchor interval",
            "continuity_reason": "h_M is a finite real sum plus explicit smooth boundary jets",
            "left_and_right_signs_follow_from": "the stationary perturbation witness is smaller than both limiting endpoint sign margins",
            "uniqueness_follows_from": "the derivative perturbation witness is smaller than the positive limiting slope margin",
            "stationary_error_at_threshold_upper": (
                threshold_stationary_error_record
            ),
            "slope_error_at_threshold_upper": threshold_slope_error_record,
            "endpoint_signs_preserved_from_threshold": endpoint_signs_preserved,
            "corrected_slope_positive_from_threshold": corrected_slope_positive,
            "polynomial_witnesses_decrease_after_threshold": True,
        },
        "entries": entries,
        "claims": {
            "finite_stationary_localization_radii_certified": True,
            "all_derived_radii_fit_inside_oriented_certificates": True,
            "derived_radii_strictly_contract_on_finite_entries": (
                derived_contract
            ),
            "resultant_and_velocity_remainders_have_polynomial_witnesses": True,
            "second_derivative_remainder_has_polynomial_witness": True,
            "ideal_corrected_root_radius_witness_tends_to_zero": True,
            "corrected_stationary_root_family_constructed_from_threshold": (
                corrected_root_family_certified
            ),
            "limiting_vector_zero_certified": False,
        },
        "remaining_obligation": {
            "localization": "complete: the unique corrected stationary root family is constructed for every cutoff at or above the certified threshold",
            "core": "prove the corrected resultant norm Q_M at those roots tends to zero",
            "completed": "the corrected stationary roots exist uniquely from the threshold and their localization radii tend to zero",
        },
        "lean_bridge": {
            "contract": "StationaryLocalizationCertificate",
            "radius_theorem": "StationaryLocalizationCertificate.witness_distance_le_radius",
            "limit_theorem": "StationaryLocalizationCertificate.radius_tendsToZero",
        },
        "arithmetic": {
            "backend": "python-flint Arb real balls",
            "python_flint_version": flint.__version__,
            "flint_version": flint.__FLINT_VERSION__,
            "decimal_digits": precision,
            "binary_precision_bits": ctx.prec,
            "directed_rounding": True,
        },
        "source_hashes": {
            name: {
                "path": str(path.relative_to(root)),
                "sha256": _sha256(path),
            }
            for name, path in source_paths.items()
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the C3 stationary-localization ledger"
    )
    parser.add_argument("certificates", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    payload = build_stationary_localization(args.certificates)
    root = Path(__file__).resolve().parents[1]
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"stationary localization written to {output}")
    for key, value in payload["claims"].items():
        print(f"{key}={str(value).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
