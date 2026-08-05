#!/usr/bin/env python3
"""Directed-rounding certification in the native real plane.

The finite operator is evaluated as

    R(t) = (X(t), Y(t)),
    F(t) = X(t)^2 + Y(t)^2,
    H(t) = F'(t) / 2 = X(t) X'(t) + Y(t) Y'(t).

For a closed interval I=[a,b], the conditions

    H(a) < 0,  H(b) > 0,  H'(I) > 0

prove that I contains a unique stationary point of F and that this point is a
strict minimum relative to I.  Coordinate enclosures for R(I) are checked
separately: excluding zero from either coordinate proves that the interval
contains no finite vector zero.

Only real Arb balls and the native R2 rotation are used.  No coordinate or
comparison theory outside the real operator is part of this module.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Iterable, Sequence

import flint
from flint import arb, ctx


SCHEMA = "org.native-carry.real-interval-minimum/v1"
OBJECTIVE = "raw_resultant_quadratic_energy"


@dataclass(frozen=True)
class SparseGeometry:
    camera: int
    cutoff: int
    geometry: str
    half_range: int
    coordinate_count: int
    largest_center: int
    terms: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class PreparedTerm:
    amplitude: Any
    log_n: Any


@dataclass(frozen=True)
class RealEvaluation:
    resultant_x: Any
    resultant_y: Any
    derivative_x: Any
    derivative_y: Any
    second_x: Any | None
    second_y: Any | None
    stationary: Any
    stationary_derivative: Any | None


def _add_coefficient(coefficients: dict[int, int], n: int, value: int) -> None:
    if n < 1:
        raise ValueError("the finite camera produced a nonpositive state index")
    coefficients[n] = coefficients.get(n, 0) + value


def build_sparse_geometry(camera: int, cutoff: int) -> SparseGeometry:
    """Build the exact integer coefficient map of the finite real operator."""
    camera = int(camera)
    cutoff = int(cutoff)
    if camera < 2:
        raise ValueError("camera must be at least 2")
    if cutoff < 1:
        raise ValueError("cutoff must be at least 1")

    coefficients: dict[int, int] = {}
    if camera == 2:
        half_range = 1
        geometry = "c2_aligned_centers_4m"
        seeds = (1,)
        centers_and_radii: Iterable[tuple[int, int]] = (
            (4 * m, 1) for m in range(1, cutoff + 1)
        )
        coordinate_count = cutoff + 1
        largest_center = 4 * cutoff
    else:
        half_range = camera // 2
        geometry = (
            "natural_saturated_even_antipode"
            if camera % 2 == 0
            else "natural_saturated_odd_width"
        )
        seeds = tuple(range(1, half_range + 1))
        centers_and_radii = (
            (camera * m, radius)
            for m in range(1, cutoff + 1)
            for radius in range(1, half_range + 1)
        )
        coordinate_count = half_range * (cutoff + 1)
        largest_center = camera * cutoff

    for seed in seeds:
        _add_coefficient(coefficients, seed, 1)
    for center, radius in centers_and_radii:
        _add_coefficient(coefficients, center - radius, 1)
        _add_coefficient(coefficients, center, -2)
        _add_coefficient(coefficients, center + radius, 1)

    terms = tuple(
        (n, coefficient)
        for n, coefficient in sorted(coefficients.items())
        if coefficient
    )
    return SparseGeometry(
        camera=camera,
        cutoff=cutoff,
        geometry=geometry,
        half_range=half_range,
        coordinate_count=coordinate_count,
        largest_center=largest_center,
        terms=terms,
    )


def prepare_terms(geometry: SparseGeometry) -> tuple[PreparedTerm, ...]:
    prepared: list[PreparedTerm] = []
    for n, coefficient in geometry.terms:
        n_ball = arb(n)
        prepared.append(
            PreparedTerm(
                amplitude=arb(coefficient) / n_ball.sqrt(),
                log_n=n_ball.log(),
            )
        )
    return tuple(prepared)


def evaluate_real_operator(
    time: Any,
    prepared: Sequence[PreparedTerm],
    *,
    second: bool,
) -> RealEvaluation:
    """Enclose R, R', H and optionally R'' and H' at a real ball."""
    rx = arb(0)
    ry = arb(0)
    dx = arb(0)
    dy = arb(0)
    ddx = arb(0) if second else None
    ddy = arb(0) if second else None

    for term in prepared:
        angle = -time * term.log_n
        cosine = angle.cos()
        sine = angle.sin()
        x = term.amplitude * cosine
        y = term.amplitude * sine
        rx += x
        ry += y
        dx += term.amplitude * term.log_n * sine
        dy -= term.amplitude * term.log_n * cosine
        if second:
            assert ddx is not None and ddy is not None
            log_sq = term.log_n * term.log_n
            ddx -= term.amplitude * log_sq * cosine
            ddy -= term.amplitude * log_sq * sine

    stationary = rx * dx + ry * dy
    stationary_derivative = None
    if second:
        assert ddx is not None and ddy is not None
        stationary_derivative = dx * dx + dy * dy + rx * ddx + ry * ddy

    return RealEvaluation(
        resultant_x=rx,
        resultant_y=ry,
        derivative_x=dx,
        derivative_y=dy,
        second_x=ddx,
        second_y=ddy,
        stationary=stationary,
        stationary_derivative=stationary_derivative,
    )


def _decimal_endpoints(center: str, radius: str) -> tuple[str, str]:
    center_value = Decimal(center)
    radius_value = Decimal(radius)
    if not center_value.is_finite() or not radius_value.is_finite():
        raise ValueError("center and radius must be finite decimals")
    if radius_value <= 0:
        raise ValueError("radius must be positive")
    digit_budget = max(len(center), len(radius)) + 32
    with localcontext() as decimal_context:
        decimal_context.prec = digit_budget
        lower = center_value - radius_value
        upper = center_value + radius_value
    return format(lower, "f"), format(upper, "f")


def _ball_payload(value: Any, digits: int) -> dict[str, Any]:
    midpoint, radius, exponent = value.mid_rad_10exp(digits)
    return {
        "midpoint_integer": str(midpoint),
        "radius_integer": str(radius),
        "exponent10": int(exponent),
        "exact_enclosure_law": "(midpoint_integer +/- radius_integer) * 10^exponent10",
        "display": value.str(digits, radius=True, more=True),
        "contains_zero": bool(value.contains(0)),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def certify_finite_minimum(
    *,
    camera: int,
    cutoff: int,
    center: str,
    radius: str,
    dps: int = 120,
    source_root: Path | None = None,
) -> dict[str, Any]:
    """Return a machine-checkable real interval certificate ledger."""
    if dps < 30:
        raise ValueError("dps must be at least 30")
    lower_text, upper_text = _decimal_endpoints(center, radius)

    ctx.dps = int(dps)
    geometry = build_sparse_geometry(camera, cutoff)
    prepared = prepare_terms(geometry)

    lower = arb(lower_text)
    upper = arb(upper_text)
    domain = arb(arb(center), arb(radius))
    lower_evaluation = evaluate_real_operator(lower, prepared, second=False)
    upper_evaluation = evaluate_real_operator(upper, prepared, second=False)
    interval_evaluation = evaluate_real_operator(domain, prepared, second=True)
    assert interval_evaluation.stationary_derivative is not None

    left_negative = bool(lower_evaluation.stationary < 0)
    right_positive = bool(upper_evaluation.stationary > 0)
    derivative_positive = bool(interval_evaluation.stationary_derivative > 0)
    unique_stationary = left_negative and right_positive and derivative_positive
    x_excludes_zero = not interval_evaluation.resultant_x.contains(0)
    y_excludes_zero = not interval_evaluation.resultant_y.contains(0)
    vector_zero_absent = x_excludes_zero or y_excludes_zero

    if not unique_stationary:
        raise RuntimeError(
            "the requested interval did not certify a unique strict stationary minimum"
        )

    root = source_root or Path(__file__).resolve().parents[1]
    source_paths = {
        "certifier": Path(__file__).resolve(),
        "finite_operator": root
        / "laboratory/native_carry_primitive_real_operator_all_bases.py",
        "precision_ladder": root / "laboratory/native_carry_precision_ladder.py",
    }
    source_hashes = {
        name: {"path": str(path.relative_to(root)), "sha256": _sha256(path)}
        for name, path in source_paths.items()
    }

    digits = dps + 12
    return {
        "schema": SCHEMA,
        "status": "REAL_INTERVAL_CERTIFIED_FINITE_MINIMUM",
        "scope": {
            "coordinate_field": "R^2",
            "objective": OBJECTIVE,
            "stationary_equation": "H = X*X' + Y*Y' = 0",
            "stationary_slope": "H' = X'^2 + Y'^2 + X*X'' + Y*Y''",
            "theory_scope": "finite_native_real_operator_only",
        },
        "operator": {
            "camera": geometry.camera,
            "cutoff": geometry.cutoff,
            "geometry": geometry.geometry,
            "half_range": geometry.half_range,
            "coordinate_count": geometry.coordinate_count,
            "largest_center": geometry.largest_center,
            "sparse_term_count": len(geometry.terms),
        },
        "domain": {
            "requested_center": center,
            "requested_radius": radius,
            "exact_decimal_lower": lower_text,
            "exact_decimal_upper": upper_text,
            "arb_enclosure": _ball_payload(domain, digits),
        },
        "enclosures": {
            "H_at_lower": _ball_payload(lower_evaluation.stationary, digits),
            "H_at_upper": _ball_payload(upper_evaluation.stationary, digits),
            "H_prime_on_domain": _ball_payload(
                interval_evaluation.stationary_derivative, digits
            ),
            "resultant_x_on_domain": _ball_payload(
                interval_evaluation.resultant_x, digits
            ),
            "resultant_y_on_domain": _ball_payload(
                interval_evaluation.resultant_y, digits
            ),
        },
        "verified_conditions": {
            "H_lower_strictly_negative": left_negative,
            "H_upper_strictly_positive": right_positive,
            "H_prime_strictly_positive_on_domain": derivative_positive,
            "resultant_x_excludes_zero": x_excludes_zero,
            "resultant_y_excludes_zero": y_excludes_zero,
        },
        "claims": {
            "unique_stationary_point_in_domain": unique_stationary,
            "strict_finite_minimum_in_domain": unique_stationary,
            "finite_vector_zero_absent_from_domain": vector_zero_absent,
            "limiting_zero_certified": False,
        },
        "logic": {
            "existence": "continuity of finite real sums plus H(lower)<0<H(upper)",
            "uniqueness": "H' is strictly positive on the full closed interval",
            "minimum": "F'=2H changes from negative to positive exactly once",
            "zero_exclusion": "at least one resultant coordinate excludes zero on the full interval",
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
        description="Certify a finite native-carry minimum with real Arb balls"
    )
    parser.add_argument("--camera", type=int, required=True)
    parser.add_argument("--cutoff", type=int, required=True)
    parser.add_argument("--center", required=True)
    parser.add_argument("--radius", required=True)
    parser.add_argument("--dps", type=int, default=120)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    certificate = certify_finite_minimum(
        camera=args.camera,
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
