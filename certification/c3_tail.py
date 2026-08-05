#!/usr/bin/env python3
"""Rigorous real tail bounds and a limiting-minimum bridge for camera C3.

For C3, the part omitted after ``cutoff = M`` is

    T_M(t) = sum_{m=M+1}^infinity Delta_1^2 psi_t(3m).

The centered second difference is written as an integral against the positive
tent kernel.  Real-plane derivative norms then give explicit majorants for
``T_M``, ``dT_M/dt`` and ``d^2 T_M/dt^2``.  The infinite sums are bounded by
closed-form real integrals of ``x^(-5/2) log(x)^q``, for ``q = 0, 1, 2``.

The same ledger combines those tail bounds with subdivided Arb evaluation of
the finite stationary numerator.  It can certify a unique strict stationary
minimum of the limiting C3 energy.  This is deliberately not a certificate
that the limiting real-plane vector vanishes at that minimum.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Sequence

import flint
from flint import arb, ctx

from certification.real_interval import (
    _ball_payload,
    _decimal_endpoints,
    build_sparse_geometry,
    evaluate_real_operator,
    prepare_terms,
)


SCHEMA = "org.native-carry.real-tail-limit-minimum/v1"


@dataclass(frozen=True)
class C3TailMajorants:
    cutoff: int
    lower_integration_point: int
    time_abs_upper: Any
    first_space_derivative_factor: Any
    second_space_derivative_factor: Any
    integral_moment_0: Any
    integral_moment_1: Any
    integral_moment_2: Any
    resultant: Any
    first_time_derivative: Any
    second_time_derivative: Any


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _positive_upper_record(value: Any, digits: int) -> dict[str, Any]:
    enclosure = _ball_payload(value, digits)
    midpoint = int(enclosure["midpoint_integer"])
    radius = int(enclosure["radius_integer"])
    upper = midpoint + radius
    if upper <= 0:
        raise RuntimeError("expected a strictly positive upper bound")
    return {
        "formula_enclosure": enclosure,
        "certified_upper": {
            "significand_integer": str(upper),
            "exponent10": int(enclosure["exponent10"]),
            "exact_value_law": "significand_integer * 10^exponent10",
        },
    }


def _upper_ball(record: dict[str, Any]) -> Any:
    upper = record["certified_upper"]
    return arb(f"{upper['significand_integer']}e{upper['exponent10']}")


def _plane_norm_upper(x_value: Any, y_value: Any, digits: int) -> dict[str, Any]:
    # abs_upper() is an outward-rounded scalar cap even when a coordinate ball
    # crosses zero.  Serializing the upper endpoint once more produces an exact
    # decimal rational used by every subsequent perturbation inequality.
    norm_cap = (
        x_value.abs_upper() * x_value.abs_upper()
        + y_value.abs_upper() * y_value.abs_upper()
    ).sqrt()
    return _positive_upper_record(norm_cap, digits)


def c3_tail_majorants(*, cutoff: int, time_abs_upper: Any) -> C3TailMajorants:
    """Return Arb enclosures of the explicit C3 tail majorants.

    The integral comparison used here requires ``cutoff >= 2``.  Then
    ``A = 3*cutoff - 1 >= 5`` and every function
    ``x^(-5/2) log(x)^q``, ``q=0,1,2``, is decreasing on ``[A, infinity)``.
    """
    cutoff = int(cutoff)
    if cutoff < 2:
        raise ValueError("the monotone C3 tail majorant requires cutoff >= 2")

    time_cap = arb(time_abs_upper)
    if time_cap < 0:
        raise ValueError("time_abs_upper must be nonnegative")

    lower_point = 3 * cutoff - 1
    lower = arb(lower_point)
    logarithm = lower.log()

    a1 = (time_cap * time_cap + arb(1) / 4).sqrt()
    a2 = (
        (time_cap * time_cap + arb(1) / 4)
        * (time_cap * time_cap + arb(9) / 4)
    ).sqrt()

    radial = lower ** (arb(-3) / 2)
    integral_0 = radial * arb(2) / 3
    integral_1 = radial * (arb(2) * logarithm / 3 + arb(4) / 9)
    integral_2 = radial * (
        arb(2) * logarithm * logarithm / 3
        + arb(8) * logarithm / 9
        + arb(16) / 27
    )

    tail_0 = a2 * integral_0 / 3
    tail_1 = (a2 * integral_1 + (2 * a1 + 1) * integral_0) / 3
    tail_2 = (
        a2 * integral_2
        + (4 * a1 + 2) * integral_1
        + 2 * integral_0
    ) / 3

    return C3TailMajorants(
        cutoff=cutoff,
        lower_integration_point=lower_point,
        time_abs_upper=time_cap,
        first_space_derivative_factor=a1,
        second_space_derivative_factor=a2,
        integral_moment_0=integral_0,
        integral_moment_1=integral_1,
        integral_moment_2=integral_2,
        resultant=tail_0,
        first_time_derivative=tail_1,
        second_time_derivative=tail_2,
    )


def _decimal_subintervals(
    lower_text: str, upper_text: str, count: int
) -> list[tuple[str, str, str, str]]:
    if count < 1:
        raise ValueError("subinterval count must be positive")
    lower = Decimal(lower_text)
    upper = Decimal(upper_text)
    with localcontext() as decimal_context:
        decimal_context.prec = max(len(lower_text), len(upper_text)) + 40
        step = (upper - lower) / Decimal(count)
        result: list[tuple[str, str, str, str]] = []
        for index in range(count):
            left = lower + Decimal(index) * step
            right = lower + Decimal(index + 1) * step
            midpoint = (left + right) / 2
            radius = (right - left) / 2
            result.append(
                tuple(format(value, "f") for value in (left, right, midpoint, radius))
            )
    return result


def certify_c3_tail_and_limit_minimum(
    *,
    cutoff: int,
    center: str,
    radius: str,
    subdivisions: int = 10,
    dps: int = 100,
    source_root: Path | None = None,
) -> dict[str, Any]:
    """Build a directed-rounding C3 tail and limiting-minimum ledger."""
    if dps < 40:
        raise ValueError("dps must be at least 40")
    if subdivisions < 1:
        raise ValueError("subdivisions must be positive")

    lower_text, upper_text = _decimal_endpoints(center, radius)
    exact_lower = Decimal(lower_text)
    exact_upper = Decimal(upper_text)
    # Decimal.copy_abs is exact and independent of the ambient decimal context.
    time_cap_text = format(
        max(exact_lower.copy_abs(), exact_upper.copy_abs()), "f"
    )

    ctx.dps = int(dps)
    digits = dps + 12
    geometry = build_sparse_geometry(3, cutoff)
    prepared = prepare_terms(geometry)
    majorants = c3_tail_majorants(
        cutoff=cutoff, time_abs_upper=arb(time_cap_text)
    )
    log_lower_point = arb(majorants.lower_integration_point).log()
    decreasing_precondition = bool(log_lower_point > arb(4) / 5)

    tail_records = {
        "resultant": _positive_upper_record(majorants.resultant, digits),
        "first_time_derivative": _positive_upper_record(
            majorants.first_time_derivative, digits
        ),
        "second_time_derivative": _positive_upper_record(
            majorants.second_time_derivative, digits
        ),
    }
    epsilon_0 = _upper_ball(tail_records["resultant"])
    epsilon_1 = _upper_ball(tail_records["first_time_derivative"])
    epsilon_2 = _upper_ball(tail_records["second_time_derivative"])

    left_evaluation = evaluate_real_operator(arb(lower_text), prepared, second=True)
    right_evaluation = evaluate_real_operator(arb(upper_text), prepared, second=True)
    center_evaluation = evaluate_real_operator(arb(center), prepared, second=True)

    def endpoint_record(evaluation: Any, *, side: str) -> dict[str, Any]:
        norm_r = _plane_norm_upper(
            evaluation.resultant_x, evaluation.resultant_y, digits
        )
        norm_d = _plane_norm_upper(
            evaluation.derivative_x, evaluation.derivative_y, digits
        )
        perturbation = (
            _upper_ball(norm_r) * epsilon_1
            + _upper_ball(norm_d) * epsilon_0
            + epsilon_0 * epsilon_1
        )
        perturbation_record = _positive_upper_record(perturbation, digits)
        perturbation_cap = _upper_ball(perturbation_record)
        margin = (
            evaluation.stationary + perturbation_cap
            if side == "left"
            else evaluation.stationary - perturbation_cap
        )
        return {
            "finite_H_enclosure": _ball_payload(evaluation.stationary, digits),
            "finite_resultant_norm_upper": norm_r,
            "finite_first_derivative_norm_upper": norm_d,
            "H_tail_perturbation_upper": perturbation_record,
            "limiting_H_sign_margin": _ball_payload(margin, digits),
        }

    left_record = endpoint_record(left_evaluation, side="left")
    right_record = endpoint_record(right_evaluation, side="right")

    center_norm_r = _plane_norm_upper(
        center_evaluation.resultant_x, center_evaluation.resultant_y, digits
    )
    center_norm_d = _plane_norm_upper(
        center_evaluation.derivative_x, center_evaluation.derivative_y, digits
    )
    center_det = (
        center_evaluation.resultant_x * center_evaluation.derivative_y
        - center_evaluation.resultant_y * center_evaluation.derivative_x
    )
    center_det_perturbation = (
        _upper_ball(center_norm_r) * epsilon_1
        + _upper_ball(center_norm_d) * epsilon_0
        + epsilon_0 * epsilon_1
    )
    center_det_perturbation_record = _positive_upper_record(
        center_det_perturbation, digits
    )
    center_limiting_det_outer = center_det + arb(
        0, _upper_ball(center_det_perturbation_record)
    )

    cover_records: list[dict[str, Any]] = []
    for index, (cell_lower, cell_upper, midpoint, cell_radius) in enumerate(
        _decimal_subintervals(lower_text, upper_text, subdivisions)
    ):
        cell = arb(arb(midpoint), arb(cell_radius))
        evaluation = evaluate_real_operator(cell, prepared, second=True)
        assert evaluation.stationary_derivative is not None
        assert evaluation.second_x is not None and evaluation.second_y is not None

        norm_r = _plane_norm_upper(
            evaluation.resultant_x, evaluation.resultant_y, digits
        )
        norm_d = _plane_norm_upper(
            evaluation.derivative_x, evaluation.derivative_y, digits
        )
        norm_dd = _plane_norm_upper(
            evaluation.second_x, evaluation.second_y, digits
        )
        perturbation = (
            2 * _upper_ball(norm_d) * epsilon_1
            + epsilon_1 * epsilon_1
            + _upper_ball(norm_r) * epsilon_2
            + epsilon_0 * _upper_ball(norm_dd)
            + epsilon_0 * epsilon_2
        )
        perturbation_record = _positive_upper_record(perturbation, digits)
        slope_margin = evaluation.stationary_derivative - _upper_ball(
            perturbation_record
        )
        cover_records.append(
            {
                "index": index,
                "exact_decimal_lower": cell_lower,
                "exact_decimal_upper": cell_upper,
                "finite_H_prime_enclosure": _ball_payload(
                    evaluation.stationary_derivative, digits
                ),
                "finite_first_derivative_x_enclosure": _ball_payload(
                    evaluation.derivative_x, digits
                ),
                "finite_resultant_norm_upper": norm_r,
                "finite_first_derivative_norm_upper": norm_d,
                "finite_second_derivative_norm_upper": norm_dd,
                "H_prime_tail_perturbation_upper": perturbation_record,
                "limiting_H_prime_lower_margin": _ball_payload(
                    slope_margin, digits
                ),
                "limiting_first_derivative_x_negative_margin": _ball_payload(
                    evaluation.derivative_x + epsilon_1, digits
                ),
            }
        )

    left_negative = (
        int(left_record["limiting_H_sign_margin"]["midpoint_integer"])
        + int(left_record["limiting_H_sign_margin"]["radius_integer"])
        < 0
    )
    right_positive = (
        int(right_record["limiting_H_sign_margin"]["midpoint_integer"])
        - int(right_record["limiting_H_sign_margin"]["radius_integer"])
        > 0
    )
    cover_positive = all(
        int(cell["limiting_H_prime_lower_margin"]["midpoint_integer"])
        - int(cell["limiting_H_prime_lower_margin"]["radius_integer"])
        > 0
        for cell in cover_records
    )
    derivative_x_negative = all(
        int(cell["limiting_first_derivative_x_negative_margin"]["midpoint_integer"])
        + int(cell["limiting_first_derivative_x_negative_margin"]["radius_integer"])
        < 0
        for cell in cover_records
    )
    limit_minimum = left_negative and right_positive and cover_positive
    if not (limit_minimum and derivative_x_negative):
        raise RuntimeError(
            "the requested domain did not certify a unique limiting minimum"
        )

    root = source_root or Path(__file__).resolve().parents[1]
    source_paths = {
        "tail_certifier": Path(__file__).resolve(),
        "finite_interval_evaluator": root / "certification/real_interval.py",
        "finite_operator": root
        / "laboratory/native_carry_primitive_real_operator_all_bases.py",
    }
    source_hashes = {
        name: {"path": str(path.relative_to(root)), "sha256": _sha256(path)}
        for name, path in source_paths.items()
    }

    return {
        "schema": SCHEMA,
        "status": "REAL_INTERVAL_CERTIFIED_C3_TAIL_AND_LIMIT_MINIMUM",
        "scope": {
            "coordinate_field": "R^2",
            "theory_scope": "native_real_operator_only",
            "limit_object": "uniform limit of the C3 centered-bracket series",
            "objective": "limiting_resultant_quadratic_energy",
            "limiting_zero_is_not_inferred_from_minimum": True,
        },
        "operator": {
            "camera": 3,
            "finite_cutoff": geometry.cutoff,
            "geometry": geometry.geometry,
            "tail_centers": "3*m for integers m >= finite_cutoff + 1",
            "tail_radius": 1,
        },
        "domain": {
            "requested_center": center,
            "requested_radius": radius,
            "exact_decimal_lower": lower_text,
            "exact_decimal_upper": upper_text,
            "time_abs_upper": time_cap_text,
            "subdivision_count": subdivisions,
        },
        "tail_derivation": {
            "bracket_identity": "Delta_1^2 f(c) = integral[-1,1] (1-|u|) * f''(c+u) du",
            "tent_kernel_mass": 1,
            "lower_integration_point": majorants.lower_integration_point,
            "log_lower_integration_point": _ball_payload(
                log_lower_point, digits
            ),
            "monotone_integral_comparison": "sum_{m=M+1}^infinity B(3m-1) <= (1/3) * integral_{3M-1}^infinity B(x) dx",
            "space_derivative_norms": {
                "first_factor_A1": _ball_payload(
                    majorants.first_space_derivative_factor, digits
                ),
                "second_factor_A2": _ball_payload(
                    majorants.second_space_derivative_factor, digits
                ),
            },
            "integral_moments": {
                "I0": _ball_payload(majorants.integral_moment_0, digits),
                "I1": _ball_payload(majorants.integral_moment_1, digits),
                "I2": _ball_payload(majorants.integral_moment_2, digits),
            },
            "verified_preconditions": {
                "cutoff_at_least_two": cutoff >= 2,
                "lower_integration_point_at_least_five": (
                    majorants.lower_integration_point >= 5
                ),
                "time_domain_is_bounded": True,
                "log_lower_point_strictly_above_four_fifths": (
                    decreasing_precondition
                ),
                "majorants_decrease_on_integration_domain": (
                    decreasing_precondition
                ),
            },
            "vanishing_with_cutoff": "each bound is (3M-1)^(-3/2) times a polynomial of degree at most two in log(3M-1)",
        },
        "tail_bounds": tail_records,
        "limit_bridge": {
            "stationary_equation": "H_infinity = R_infinity dot R_infinity'",
            "stationary_derivative": "H_infinity' = |R_infinity'|^2 + R_infinity dot R_infinity''",
            "left_endpoint": left_record,
            "right_endpoint": right_record,
            "positive_slope_cover": cover_records,
        },
        "vector_zero_reduction": {
            "oriented_determinant": "K_infinity = det(R_infinity, R_infinity') = X_infinity*Y_infinity' - Y_infinity*X_infinity'",
            "real_plane_identity": "H_infinity^2 + K_infinity^2 = |R_infinity|^2 * |R_infinity'|^2",
            "certified_nonvanishing_velocity": "the x coordinate of R_infinity' is strictly negative on every cover cell",
            "at_unique_stationary_point": "R_infinity = 0 iff K_infinity = 0",
            "finite_cutoff_reference_at_requested_center": {
                "finite_K_enclosure": _ball_payload(center_det, digits),
                "K_tail_perturbation_upper": center_det_perturbation_record,
                "limiting_K_outer_enclosure": _ball_payload(
                    center_limiting_det_outer, digits
                ),
                "current_norm_tail_bound_decides_determinant_sign": (
                    not center_limiting_det_outer.contains(0)
                ),
            },
            "determinant_zero_certified": False,
        },
        "verified_conditions": {
            "limiting_H_left_strictly_negative": left_negative,
            "limiting_H_right_strictly_positive": right_positive,
            "limiting_H_prime_strictly_positive_on_cover": cover_positive,
            "limiting_first_derivative_x_strictly_negative_on_cover": (
                derivative_x_negative
            ),
        },
        "claims": {
            "uniform_C3_tail_bound_for_resultant": True,
            "uniform_C3_tail_bound_for_first_time_derivative": True,
            "uniform_C3_tail_bound_for_second_time_derivative": True,
            "C3_limit_is_twice_continuously_differentiable_on_domain": True,
            "unique_limiting_stationary_point_in_domain": limit_minimum,
            "strict_limiting_minimum_in_domain": limit_minimum,
            "limiting_first_derivative_nonzero_on_domain": (
                derivative_x_negative
            ),
            "vector_zero_reduced_to_determinant_at_stationary_point": (
                limit_minimum and derivative_x_negative
            ),
            "limiting_vector_zero_certified": False,
        },
        "logic": {
            "regularity": "uniform tail bounds through order two justify termwise time differentiation and continuity",
            "existence": "continuity plus H_infinity(lower) < 0 < H_infinity(upper)",
            "uniqueness": "the subdivided cover proves H_infinity' > 0 throughout the domain",
            "minimum": "the unique stationary numerator changes from negative to positive",
            "zero_boundary": "stationarity of the limiting energy does not imply simultaneous vanishing of both real coordinates",
        },
        "arithmetic": {
            "backend": "python-flint Arb real balls",
            "python_flint_version": flint.__version__,
            "flint_version": flint.__FLINT_VERSION__,
            "decimal_digits": dps,
            "binary_precision_bits": ctx.prec,
            "directed_rounding": True,
        },
        "source_hashes": source_hashes,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Certify real C3 tail bounds and a limiting stationary minimum"
    )
    parser.add_argument("--cutoff", type=int, required=True)
    parser.add_argument("--center", required=True)
    parser.add_argument("--radius", required=True)
    parser.add_argument("--subdivisions", type=int, default=10)
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    certificate = certify_c3_tail_and_limit_minimum(
        cutoff=args.cutoff,
        center=args.center,
        radius=args.radius,
        subdivisions=args.subdivisions,
        dps=args.dps,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(certificate, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"certificate written to {args.output}")
    for key, value in certificate["claims"].items():
        print(f"{key}={str(value).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
