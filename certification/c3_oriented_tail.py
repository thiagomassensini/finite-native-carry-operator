#!/usr/bin/env python3
"""Oriented real C3 tail enclosures with a sixth-derivative remainder.

The norm-only tail certificate deliberately discards rotation and cancellation.
This module retains both real coordinates.  For ``F_k = d^k psi / dt^k`` and
``C = 3(M+1)``, the certified approximation is

    sum_{m>M} Delta_1^2 F_k(3m)
      = -F_k'(C)/3 + F_k''(C)/2 - 5 F_k'''(C)/18
        + F_k''''(C)/24 + F_k'''''(C)/60 + remainder.

The formula follows from the symmetric Taylor expansion of the bracket through
the fourth space derivative and real Euler--Maclaurin bounds for the two
resulting lattice sums.  The remainder uses only an integral majorant for the
sixth space derivative.  Everything is evaluated in the native real plane.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from decimal import Decimal
from math import comb, factorial
from pathlib import Path
from typing import Any, Sequence

import flint
from flint import arb, ctx

from certification.c3_tail import _plane_norm_upper, _positive_upper_record, _upper_ball
from certification.real_interval import (
    _ball_payload,
    _decimal_endpoints,
    build_sparse_geometry,
    evaluate_real_operator,
    prepare_terms,
)


SCHEMA = "org.native-carry.real-oriented-tail-limit-minimum/v1"


@dataclass(frozen=True)
class OrientedTailEnclosure:
    time_derivative_order: int
    approximation_x: Any
    approximation_y: Any
    remainder_norm: Any
    enclosure_x: Any
    enclosure_y: Any


@dataclass(frozen=True)
class LimitEvaluation:
    resultant_x: Any
    resultant_y: Any
    derivative_x: Any
    derivative_y: Any
    second_x: Any
    second_y: Any
    stationary: Any
    stationary_derivative: Any
    determinant: Any
    tails: tuple[OrientedTailEnclosure, ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _pair_add(left: tuple[Any, Any], right: tuple[Any, Any]) -> tuple[Any, Any]:
    return left[0] + right[0], left[1] + right[1]


def _pair_scale(scalar: Any, value: tuple[Any, Any]) -> tuple[Any, Any]:
    return scalar * value[0], scalar * value[1]


def _apply_j(value: tuple[Any, Any], order: int) -> tuple[Any, Any]:
    x_value, y_value = value
    for _ in range(order % 4):
        x_value, y_value = -y_value, x_value
    return x_value, y_value


def _base_space_derivative_pair(space_order: int, time: Any) -> tuple[Any, Any]:
    """Rotating-frame coefficient of the requested space derivative."""
    value = (arb(1), arb(0))
    for index in range(space_order):
        radial = -arb(1) / 2 - index
        # (radial * I - time * J) value
        value = (
            radial * value[0] + time * value[1],
            radial * value[1] - time * value[0],
        )
    return value


def _harmonic(index: int) -> Any:
    return sum((arb(1) / value for value in range(1, index + 1)), arb(0))


def _log_power_derivative(
    time_order: int, derivative_order: int, x_value: Any, logarithm: Any
) -> Any:
    """Space derivative of ``(-log x)^time_order`` for orders zero to two."""
    if time_order == 0:
        return arb(1) if derivative_order == 0 else arb(0)
    if time_order == 1:
        if derivative_order == 0:
            return -logarithm
        return (
            arb((-1) ** derivative_order * factorial(derivative_order - 1))
            * x_value ** (-derivative_order)
        )
    if time_order == 2:
        if derivative_order == 0:
            return logarithm * logarithm
        return (
            2
            * arb((-1) ** (derivative_order - 1) * factorial(derivative_order - 1))
            * (logarithm - _harmonic(derivative_order - 1))
            * x_value ** (-derivative_order)
        )
    raise ValueError("only time derivative orders zero, one and two are supported")


def mixed_space_time_derivative(
    *, x_value: Any, time: Any, space_order: int, time_order: int
) -> tuple[Any, Any]:
    """Enclose ``d_x^space_order d_t^time_order psi_t(x)`` in R2."""
    if not 0 <= time_order <= 2:
        raise ValueError("time_order must be zero, one or two")
    if space_order < 0:
        raise ValueError("space_order must be nonnegative")

    x_ball = arb(x_value)
    time_ball = arb(time)
    logarithm = x_ball.log()
    angle = -time_ball * logarithm
    rotation = (angle.cos(), angle.sin())
    result = (arb(0), arb(0))

    for derivative_order in range(space_order + 1):
        base_order = space_order - derivative_order
        frame = _base_space_derivative_pair(base_order, time_ball)
        rotated = (
            rotation[0] * frame[0] - rotation[1] * frame[1],
            rotation[1] * frame[0] + rotation[0] * frame[1],
        )
        base = _pair_scale(
            x_ball ** (-arb(1) / 2 - base_order), rotated
        )
        base = _apply_j(base, time_order)
        coefficient = comb(space_order, derivative_order) * _log_power_derivative(
            time_order, derivative_order, x_ball, logarithm
        )
        result = _pair_add(result, _pair_scale(coefficient, base))
    return result


def _space_factor(space_order: int, time_abs_upper: Any) -> Any:
    factor = arb(1)
    time_cap = arb(time_abs_upper)
    for index in range(space_order):
        half_integer = arb(index) + arb(1) / 2
        factor *= (time_cap * time_cap + half_integer * half_integer).sqrt()
    return factor


def _sixth_integral_moments(lower: Any) -> tuple[Any, Any, Any]:
    lower_ball = arb(lower)
    logarithm = lower_ball.log()
    radial = lower_ball ** (-arb(11) / 2)
    return (
        radial * arb(2) / 11,
        radial * (arb(2) * logarithm / 11 + arb(4) / 121),
        radial
        * (
            arb(2) * logarithm * logarithm / 11
            + arb(8) * logarithm / 121
            + arb(16) / 1331
        ),
    )


def sixth_derivative_integral_majorant(
    *, lower: Any, time_abs_upper: Any, time_order: int
) -> Any:
    """Bound the integral of ``|d_x^6 d_t^k psi|`` from lower to infinity."""
    if not 0 <= time_order <= 2:
        raise ValueError("time_order must be zero, one or two")
    factors = tuple(
        _space_factor(order, time_abs_upper) for order in range(7)
    )
    moment_0, moment_1, moment_2 = _sixth_integral_moments(lower)
    if time_order == 0:
        return factors[6] * moment_0
    if time_order == 1:
        constant = (
            6 * factors[5]
            + 15 * factors[4]
            + 40 * factors[3]
            + 90 * factors[2]
            + 144 * factors[1]
            + 120
        )
        return factors[6] * moment_1 + constant * moment_0

    log_coefficient = (
        12 * factors[5]
        + 30 * factors[4]
        + 80 * factors[3]
        + 180 * factors[2]
        + 288 * factors[1]
        + 240
    )
    constant = (
        30 * factors[4]
        + 120 * factors[3]
        + 330 * factors[2]
        + 600 * factors[1]
        + 548
    )
    return (
        factors[6] * moment_2
        + log_coefficient * moment_1
        + constant * moment_0
    )


def oriented_c3_tail_enclosure(
    *, cutoff: int, time: Any, time_abs_upper: Any, time_order: int
) -> OrientedTailEnclosure:
    """Enclose one oriented C3 tail derivative by its fifth-order boundary jet."""
    cutoff = int(cutoff)
    if cutoff < 2:
        raise ValueError("the oriented C3 tail requires cutoff >= 2")
    if not 0 <= time_order <= 2:
        raise ValueError("time_order must be zero, one or two")

    center = arb(3 * (cutoff + 1))
    coefficients = (
        -arb(1) / 3,
        arb(1) / 2,
        -arb(5) / 18,
        arb(1) / 24,
        arb(1) / 60,
    )
    approximation = (arb(0), arb(0))
    for space_order, coefficient in enumerate(coefficients, start=1):
        derivative = mixed_space_time_derivative(
            x_value=center,
            time=time,
            space_order=space_order,
            time_order=time_order,
        )
        approximation = _pair_add(
            approximation, _pair_scale(coefficient, derivative)
        )

    euler_lower = center
    taylor_lower = arb(3 * cutoff - 1)
    remainder = (
        arb(7)
        / 120
        * sixth_derivative_integral_majorant(
            lower=euler_lower,
            time_abs_upper=time_abs_upper,
            time_order=time_order,
        )
        + arb(1)
        / 1080
        * sixth_derivative_integral_majorant(
            lower=taylor_lower,
            time_abs_upper=time_abs_upper,
            time_order=time_order,
        )
    )
    return OrientedTailEnclosure(
        time_derivative_order=time_order,
        approximation_x=approximation[0],
        approximation_y=approximation[1],
        remainder_norm=remainder,
        enclosure_x=arb(approximation[0], remainder),
        enclosure_y=arb(approximation[1], remainder),
    )


def evaluate_oriented_c3_limit(
    *, time: Any, prepared: Sequence[Any], cutoff: int, time_abs_upper: Any
) -> LimitEvaluation:
    finite = evaluate_real_operator(time, prepared, second=True)
    assert finite.second_x is not None and finite.second_y is not None
    tails = tuple(
        oriented_c3_tail_enclosure(
            cutoff=cutoff,
            time=time,
            time_abs_upper=time_abs_upper,
            time_order=order,
        )
        for order in range(3)
    )
    resultant_x = finite.resultant_x + tails[0].enclosure_x
    resultant_y = finite.resultant_y + tails[0].enclosure_y
    derivative_x = finite.derivative_x + tails[1].enclosure_x
    derivative_y = finite.derivative_y + tails[1].enclosure_y
    second_x = finite.second_x + tails[2].enclosure_x
    second_y = finite.second_y + tails[2].enclosure_y
    stationary = resultant_x * derivative_x + resultant_y * derivative_y
    stationary_derivative = (
        derivative_x * derivative_x
        + derivative_y * derivative_y
        + resultant_x * second_x
        + resultant_y * second_y
    )
    determinant = resultant_x * derivative_y - resultant_y * derivative_x
    return LimitEvaluation(
        resultant_x=resultant_x,
        resultant_y=resultant_y,
        derivative_x=derivative_x,
        derivative_y=derivative_y,
        second_x=second_x,
        second_y=second_y,
        stationary=stationary,
        stationary_derivative=stationary_derivative,
        determinant=determinant,
        tails=tails,
    )


def _tail_payload(tail: OrientedTailEnclosure, digits: int) -> dict[str, Any]:
    return {
        "time_derivative_order": tail.time_derivative_order,
        "approximation_x": _ball_payload(tail.approximation_x, digits),
        "approximation_y": _ball_payload(tail.approximation_y, digits),
        "remainder_norm_upper": _positive_upper_record(
            tail.remainder_norm, digits
        ),
        "oriented_enclosure_x": _ball_payload(tail.enclosure_x, digits),
        "oriented_enclosure_y": _ball_payload(tail.enclosure_y, digits),
    }


def certify_oriented_c3_limit_minimum(
    *,
    cutoff: int,
    center: str,
    radius: str,
    dps: int = 110,
    source_root: Path | None = None,
) -> dict[str, Any]:
    """Build an oriented tail certificate around the limiting C3 minimum."""
    if dps < 50:
        raise ValueError("dps must be at least 50")
    lower_text, upper_text = _decimal_endpoints(center, radius)
    exact_lower = Decimal(lower_text)
    exact_upper = Decimal(upper_text)
    time_cap_text = format(
        max(exact_lower.copy_abs(), exact_upper.copy_abs()), "f"
    )

    ctx.dps = int(dps)
    digits = dps + 12
    geometry = build_sparse_geometry(3, cutoff)
    prepared = prepare_terms(geometry)
    time_cap = arb(time_cap_text)

    lower_evaluation = evaluate_oriented_c3_limit(
        time=arb(lower_text),
        prepared=prepared,
        cutoff=cutoff,
        time_abs_upper=time_cap,
    )
    upper_evaluation = evaluate_oriented_c3_limit(
        time=arb(upper_text),
        prepared=prepared,
        cutoff=cutoff,
        time_abs_upper=time_cap,
    )
    domain = arb(arb(center), arb(radius))
    domain_evaluation = evaluate_oriented_c3_limit(
        time=domain,
        prepared=prepared,
        cutoff=cutoff,
        time_abs_upper=time_cap,
    )
    center_evaluation = evaluate_oriented_c3_limit(
        time=arb(center),
        prepared=prepared,
        cutoff=cutoff,
        time_abs_upper=time_cap,
    )

    left_negative = bool(lower_evaluation.stationary < 0)
    right_positive = bool(upper_evaluation.stationary > 0)
    slope_positive = bool(domain_evaluation.stationary_derivative > 0)
    derivative_x_negative = bool(domain_evaluation.derivative_x < 0)
    unique_minimum = left_negative and right_positive and slope_positive
    if not (unique_minimum and derivative_x_negative):
        raise RuntimeError(
            "the requested oriented interval did not certify the limit minimum"
        )

    center_norm = _plane_norm_upper(
        center_evaluation.resultant_x, center_evaluation.resultant_y, digits
    )
    derivative_norm = _plane_norm_upper(
        domain_evaluation.derivative_x, domain_evaluation.derivative_y, digits
    )
    stationary_resultant_bound = (
        _upper_ball(center_norm)
        + _upper_ball(derivative_norm) * arb(radius)
    )
    stationary_resultant_record = _positive_upper_record(
        stationary_resultant_bound, digits
    )
    stationary_determinant_bound = (
        _upper_ball(stationary_resultant_record) * _upper_ball(derivative_norm)
    )
    stationary_determinant_record = _positive_upper_record(
        stationary_determinant_bound, digits
    )

    root = source_root or Path(__file__).resolve().parents[1]
    source_paths = {
        "oriented_tail_certifier": Path(__file__).resolve(),
        "finite_interval_evaluator": root / "certification/real_interval.py",
        "norm_tail_certifier": root / "certification/c3_tail.py",
        "finite_operator": root
        / "laboratory/native_carry_primitive_real_operator_all_bases.py",
    }
    source_hashes = {
        name: {"path": str(path.relative_to(root)), "sha256": _sha256(path)}
        for name, path in source_paths.items()
    }

    return {
        "schema": SCHEMA,
        "status": "REAL_INTERVAL_CERTIFIED_C3_ORIENTED_TAIL_LIMIT_MINIMUM",
        "scope": {
            "coordinate_field": "R^2",
            "theory_scope": "native_real_operator_only",
            "limit_object": "C3 seed plus all radius-one centered brackets at centers 3*m for integers m >= 1",
            "objective": "limiting_resultant_quadratic_energy",
            "oriented_tail": True,
        },
        "operator": {
            "camera": 3,
            "finite_cutoff": geometry.cutoff,
            "geometry": geometry.geometry,
            "tail_radius": 1,
            "first_omitted_center": 3 * (cutoff + 1),
        },
        "domain": {
            "requested_center": center,
            "requested_radius": radius,
            "exact_decimal_lower": lower_text,
            "exact_decimal_upper": upper_text,
            "time_abs_upper": time_cap_text,
            "arb_enclosure": _ball_payload(domain, digits),
        },
        "oriented_tail_theorem": {
            "boundary_jet": "-F_k'(C)/3 + F_k''(C)/2 - 5*F_k'''(C)/18 + F_k''''(C)/24 + F_k'''''(C)/60",
            "boundary_center": 3 * (cutoff + 1),
            "bracket_expansion": "Delta_1^2 F_k(c) = F_k''(c) + F_k''''(c)/12 + rho_k(c)",
            "bracket_remainder": "|rho_k(c)| <= sup_[c-1,c+1] |F_k^(6)| / 360",
            "combined_remainder": "(7/120)*integral_C^infinity |F_k^(6)| + (1/1080)*integral_(3M-1)^infinity |F_k^(6)|",
            "verified_preconditions": {
                "cutoff_at_least_two": cutoff >= 2,
                "sixth_derivative_majorants_integrable": True,
                "majorants_decrease_on_tail_domain": True,
            },
        },
        "tail_enclosures_on_domain": {
            "resultant": _tail_payload(domain_evaluation.tails[0], digits),
            "first_time_derivative": _tail_payload(
                domain_evaluation.tails[1], digits
            ),
            "second_time_derivative": _tail_payload(
                domain_evaluation.tails[2], digits
            ),
        },
        "limit_enclosures": {
            "H_at_lower": _ball_payload(lower_evaluation.stationary, digits),
            "H_at_upper": _ball_payload(upper_evaluation.stationary, digits),
            "H_prime_on_domain": _ball_payload(
                domain_evaluation.stationary_derivative, digits
            ),
            "first_derivative_x_on_domain": _ball_payload(
                domain_evaluation.derivative_x, digits
            ),
            "resultant_x_at_center": _ball_payload(
                center_evaluation.resultant_x, digits
            ),
            "resultant_y_at_center": _ball_payload(
                center_evaluation.resultant_y, digits
            ),
            "determinant_at_center": _ball_payload(
                center_evaluation.determinant, digits
            ),
        },
        "stationary_point_bounds": {
            "resultant_norm_upper": stationary_resultant_record,
            "determinant_abs_upper": stationary_determinant_record,
            "derivation": "mean-value bound from the center to the unique stationary point, followed by |det(R,R')| <= |R|*|R'|",
        },
        "verified_conditions": {
            "limiting_H_left_strictly_negative": left_negative,
            "limiting_H_right_strictly_positive": right_positive,
            "limiting_H_prime_strictly_positive_on_domain": slope_positive,
            "limiting_first_derivative_x_strictly_negative_on_domain": (
                derivative_x_negative
            ),
        },
        "claims": {
            "oriented_C3_tail_enclosed_through_second_time_derivative": True,
            "unique_limiting_stationary_point_in_domain": unique_minimum,
            "strict_limiting_minimum_in_domain": unique_minimum,
            "limiting_first_derivative_nonzero_on_domain": (
                derivative_x_negative
            ),
            "stationary_resultant_has_certified_small_norm": True,
            "stationary_determinant_has_certified_small_absolute_value": True,
            "limiting_vector_zero_certified": False,
        },
        "logic": {
            "existence": "H_infinity changes from negative to positive",
            "uniqueness": "H_infinity' is positive on the complete interval",
            "small_resultant": "the center enclosure plus the derivative Lipschitz bound controls R_infinity at the stationary point",
            "zero_boundary": "a positive upper bound, however small, is not an equality proof",
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
        description="Certify the oriented real C3 tail near its limiting minimum"
    )
    parser.add_argument("--cutoff", type=int, required=True)
    parser.add_argument("--center", required=True)
    parser.add_argument("--radius", required=True)
    parser.add_argument("--dps", type=int, default=110)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    certificate = certify_oriented_c3_limit_minimum(
        cutoff=args.cutoff,
        center=args.center,
        radius=args.radius,
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
