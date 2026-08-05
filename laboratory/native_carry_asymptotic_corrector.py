#!/usr/bin/env python3
"""Native asymptotic tail corrector for finite-cutoff carry minima.

This program does not evaluate or modify the native-carry operator.  It reads
finite-cutoff minima already emitted by ``native_carry_precision_ladder.py``
and removes the leading dyadic tail modes predicted by the centered second
-difference geometry of the native amplitude n^(-1/2).

Native tail model
-----------------
Write the state using only its native real rotation,

    f_t(n) = n^(-1/2) Rot(-t log n) e1.

A centered bracket has the expansion

    f(n-r) - 2 f(n) + f(n+r)
      = r^2 f''(n) + r^4 f''''(n)/12 + ... .

After summation beyond a cutoff M, the real vector tail has layers

    M^(-(3/2+j)) Rot(-t log M) v_j(t),    j = 0, 1, 2, ... .

Under M -> 2M, layer j is contracted and rotated by the real-plane map

    q_j(t) = 2^(-(3/2+j)) Rot(-t log 2).

For a non-degenerate limiting stationary minimum, the finite-root drift
inherits these leading damped rotational modes.  Each real mode is annihilated
by

    P_j(L) = L^2 - 2 rho_j cos(T log 2) L + rho_j^2,
    rho_j = 2^(-(3/2+j)),

where T is the limiting center and L advances one dyadic cutoff.  Multiplying
the first m factors gives a recurrence of order 2m.  A phase-consistent fixed
point applied to 2m+1 consecutive minima estimates T without Newton and without
interval arithmetic.

Canonical use
-------------
Two real rotational layers (j=0,1), with radial decays M^(-3/2) and
M^(-5/2), are the canonical corrected model.

It needs five consecutive dyadic cutoffs and can make a genuine rolling
one-step holdout prediction for the sixth.  More layers are printed only as an
exploratory linear-tail cross-check: nonlinear root-shift terms can enter before
higher linear layers, so they are not promoted automatically.

Methodological status
---------------------
* the finite minima remain the original operator outputs;
* no Newton iteration is used;
* no interval arithmetic is used;
* no rigorous remainder bound is claimed;
* rolling holdouts use only earlier cutoffs;
* the reported uncertainty is empirical model spread, not certification.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import statistics
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

try:
    import mpmath as mp
except Exception as exc:  # pragma: no cover - dependency path
    raise SystemExit(
        "mpmath is required. Install it with: python3 -m pip install mpmath"
    ) from exc


SCHEMA = "org.native-carry.asymptotic-tail-corrector/v1"
DEFAULT_DECAY = "1.5"


@dataclass(frozen=True)
class CutoffPoint:
    cutoff: int
    root_text: str
    score_text: str | None
    resultant_norm_text: str | None
    source: str
    root_digits: int


@dataclass(frozen=True)
class CenterSolve:
    center: Any
    coefficients: tuple[Any, ...]
    fixed_point_iterations: int
    fixed_point_residual: Any
    recurrence_residual: Any


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mp_text(value: Any, digits: int) -> str:
    return mp.nstr(mp.mpf(value), n=max(2, digits), strip_zeros=True)


def scientific_text(value: Any, digits: int = 18) -> str:
    value = mp.mpf(value)
    if value == 0:
        return "0"
    return mp.nstr(value, n=max(2, digits), min_fixed=0, max_fixed=0)


def count_numeric_digits(text: str) -> int:
    return sum(character.isdigit() for character in text.lstrip("+-"))


def is_power_of_two(value: int) -> bool:
    return value > 0 and (value & (value - 1)) == 0


def power_of_two_exponent(value: int) -> int:
    if not is_power_of_two(value):
        raise ValueError(f"{value} is not a power of two")
    return value.bit_length() - 1


def fixed_decimal(value: Any, places: int) -> str:
    """Round an mp number to a fixed number of decimal places without float."""
    value = mp.mpf(value)
    places = max(0, int(places))
    scale = mp.power(10, places)
    scaled = int(mp.nint(value * scale))
    sign = "-" if scaled < 0 else ""
    digits = str(abs(scaled))
    if places == 0:
        return sign + digits
    digits = digits.rjust(places + 1, "0")
    return sign + digits[:-places] + "." + digits[-places:]


def safe_decimal_places(radius: Any, cap: int = 80) -> int | None:
    """Decimal places whose half-unit exceeds the empirical radius."""
    radius = abs(mp.mpf(radius))
    if radius == 0:
        return cap
    if radius >= mp.mpf("0.5"):
        return 0
    places = int(mp.floor(-mp.log10(2 * radius)))
    return max(0, min(cap, places))


def median_mpf(values: Sequence[Any]) -> Any:
    ordered = sorted(mp.mpf(value) for value in values)
    if not ordered:
        raise ValueError("median requires at least one value")
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def convolve(first: Sequence[Any], second: Sequence[Any]) -> list[Any]:
    result = [mp.mpf("0")] * (len(first) + len(second) - 1)
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            result[i + j] += left * right
    return result


def annihilator_coefficients(
    center: Any,
    decay: Any,
    linear_layers: int,
) -> tuple[Any, ...]:
    """Return descending coefficients of product_j P_j(L)."""
    center = mp.mpf(center)
    decay = mp.mpf(decay)
    if linear_layers < 1:
        raise ValueError("linear_layers must be positive")
    cosine = mp.cos(center * mp.log(2))
    polynomial: list[Any] = [mp.mpf("1")]
    for layer in range(linear_layers):
        rho = mp.power(2, -(decay + layer))
        factor = [mp.mpf("1"), -2 * rho * cosine, rho * rho]
        polynomial = convolve(polynomial, factor)
    return tuple(polynomial)


def recurrence_residual(
    values: Sequence[Any],
    center: Any,
    coefficients: Sequence[Any],
) -> Any:
    order = len(coefficients) - 1
    if len(values) != order + 1:
        raise ValueError("wrong number of values for recurrence residual")
    errors = [mp.mpf(value) - center for value in values]
    return mp.fsum(
        coefficients[index] * errors[order - index]
        for index in range(order + 1)
    )


def phase_consistent_center(
    values: Sequence[Any],
    *,
    decay: Any,
    linear_layers: int,
    tolerance: Any,
    max_iterations: int = 200,
) -> CenterSolve:
    """Solve the annihilated-tail center by a phase-consistent fixed point.

    This is not Newton.  At a provisional T, the recurrence coefficients are
    frozen, the constant sequence is eliminated in closed form, and T is
    updated.  The map is strongly contractive in the intended asymptotic regime.
    """
    order = 2 * linear_layers
    if len(values) != order + 1:
        raise ValueError(
            f"{linear_layers} layer(s) require {order + 1} consecutive values"
        )
    values = [mp.mpf(value) for value in values]
    tolerance = abs(mp.mpf(tolerance))
    center = values[-1]
    previous_step: Any | None = None

    for iteration in range(1, max_iterations + 1):
        coefficients = annihilator_coefficients(center, decay, linear_layers)
        denominator = mp.fsum(coefficients)
        if abs(denominator) <= mp.eps:
            raise RuntimeError("annihilator is nearly singular at L=1")
        numerator = mp.fsum(
            coefficients[index] * values[order - index]
            for index in range(order + 1)
        )
        candidate = numerator / denominator
        step = candidate - center

        # Safeguard only if the raw fixed-point step starts expanding.
        if previous_step is not None and abs(step) > 1.25 * abs(previous_step):
            candidate = (candidate + center) / 2
            step = candidate - center

        center = candidate
        if abs(step) <= tolerance:
            coefficients = annihilator_coefficients(center, decay, linear_layers)
            residual = recurrence_residual(values, center, coefficients)
            return CenterSolve(
                center=center,
                coefficients=coefficients,
                fixed_point_iterations=iteration,
                fixed_point_residual=step,
                recurrence_residual=residual,
            )
        previous_step = step

    raise RuntimeError(
        "phase-consistent fixed point did not converge; reduce the number of "
        "layers or use later cutoffs"
    )


def predict_next(
    history: Sequence[Any],
    *,
    center: Any,
    coefficients: Sequence[Any],
) -> Any:
    """Advance one dyadic cutoff using the annihilator recurrence."""
    order = len(coefficients) - 1
    if len(history) != order:
        raise ValueError(f"prediction requires the last {order} values")
    center = mp.mpf(center)
    errors = [mp.mpf(value) - center for value in history]
    next_error = -mp.fsum(
        coefficients[index] * errors[order - index]
        for index in range(1, order + 1)
    )
    return center + next_error


def parse_point(text: str) -> tuple[int, str]:
    try:
        cutoff_text, root_text = text.split(":", 1)
        cutoff = int(cutoff_text.strip())
        root = root_text.strip()
        if cutoff < 1 or not root:
            raise ValueError
        mp.mpf(root)
        return cutoff, root
    except Exception as exc:
        raise ValueError(
            f"invalid --point {text!r}; expected CUTOFF:ROOT"
        ) from exc


def load_points(
    paths: Sequence[Path],
    manual_points: Sequence[str],
    *,
    requested_camera: int | None,
) -> tuple[dict[int, CutoffPoint], list[dict[str, Any]], dict[str, Any]]:
    points: dict[int, CutoffPoint] = {}
    conflicts: list[dict[str, Any]] = []
    source_meta: list[dict[str, Any]] = []
    detected_cameras: set[int] = set()
    detected_geometries: set[str] = set()
    detected_states: set[str] = set()

    def consume(point: CutoffPoint) -> None:
        previous = points.get(point.cutoff)
        if previous is None:
            points[point.cutoff] = point
            return
        work_dps = max(previous.root_digits, point.root_digits, 50) + 10
        with mp.workdps(work_dps):
            difference = abs(mp.mpf(previous.root_text) - mp.mpf(point.root_text))
            comparison_digits = min(previous.root_digits, point.root_digits)
            tolerance = mp.power(10, -max(8, comparison_digits - 8))
        if difference > tolerance:
            conflicts.append(
                {
                    "cutoff": point.cutoff,
                    "first_source": previous.source,
                    "first_root": previous.root_text,
                    "second_source": point.source,
                    "second_root": point.root_text,
                    "absolute_difference": scientific_text(difference, 20),
                }
            )
        # Prefer the textual witness carrying more digits.
        if point.root_digits > previous.root_digits:
            points[point.cutoff] = point

    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            raise ValueError(f"could not read JSON {path}: {exc}") from exc

        operator = payload.get("operator") or {}
        camera = operator.get("camera")
        if camera is not None:
            camera = int(camera)
            detected_cameras.add(camera)
            if requested_camera is not None and camera != requested_camera:
                raise ValueError(
                    f"{path} has camera {camera}, expected {requested_camera}"
                )
        geometry = operator.get("geometry")
        if geometry:
            detected_geometries.add(str(geometry))
        state = operator.get("state")
        if state:
            detected_states.add(str(state))

        rows = payload.get("cutoff_ladder")
        if not isinstance(rows, list):
            raise ValueError(f"{path} has no cutoff_ladder array")
        accepted = 0
        for row in rows:
            try:
                cutoff = int(row["cutoff"])
                root_text = str(row["root"])
                mp.mpf(root_text)
            except Exception as exc:
                raise ValueError(f"invalid cutoff row in {path}: {row}") from exc
            consume(
                CutoffPoint(
                    cutoff=cutoff,
                    root_text=root_text,
                    score_text=(str(row["score"]) if row.get("score") is not None else None),
                    resultant_norm_text=(
                        str(row["resultant_norm"])
                        if row.get("resultant_norm") is not None
                        else None
                    ),
                    source=str(path),
                    root_digits=count_numeric_digits(root_text),
                )
            )
            accepted += 1
        source_meta.append(
            {
                "path": str(path),
                "schema": payload.get("schema"),
                "status": payload.get("status"),
                "rows": accepted,
                "sha256": file_sha256(path),
            }
        )

    for text in manual_points:
        cutoff, root_text = parse_point(text)
        consume(
            CutoffPoint(
                cutoff=cutoff,
                root_text=root_text,
                score_text=None,
                resultant_norm_text=None,
                source="manual --point",
                root_digits=count_numeric_digits(root_text),
            )
        )

    if not points:
        raise ValueError("no cutoff points were loaded")
    if len(detected_cameras) > 1:
        raise ValueError(f"mixed cameras are not allowed: {sorted(detected_cameras)}")
    if len(detected_geometries) > 1:
        raise ValueError(f"mixed geometries are not allowed: {sorted(detected_geometries)}")
    if len(detected_states) > 1:
        raise ValueError("mixed operator state definitions are not allowed")

    metadata = {
        "sources": source_meta,
        "detected_camera": next(iter(detected_cameras), requested_camera),
        "detected_geometry": next(iter(detected_geometries), None),
        "detected_state": next(iter(detected_states), None),
    }
    return points, conflicts, metadata


def build_layer_report(
    dyadic_roots: dict[int, Any],
    *,
    decay: Any,
    linear_layers: int,
    dps: int,
    consensus_windows: int,
    recent_holdouts: int,
    display_digits: int,
) -> dict[str, Any]:
    order = 2 * linear_layers
    exponents = sorted(dyadic_roots)
    exponent_set = set(exponents)
    tolerance = mp.power(10, -(dps - 20))
    windows: list[dict[str, Any]] = []

    if exponents:
        for start in range(exponents[0], exponents[-1] - order + 1):
            required = list(range(start, start + order + 1))
            if not all(exponent in exponent_set for exponent in required):
                continue
            values = [dyadic_roots[exponent] for exponent in required]
            solved = phase_consistent_center(
                values,
                decay=decay,
                linear_layers=linear_layers,
                tolerance=tolerance,
            )
            coefficient_strings = [mp_text(value, display_digits) for value in solved.coefficients]
            row: dict[str, Any] = {
                "start_exponent": start,
                "cutoffs": [2**exponent for exponent in required],
                "center": mp_text(solved.center, display_digits),
                "fixed_point_iterations": solved.fixed_point_iterations,
                "fixed_point_residual": scientific_text(
                    solved.fixed_point_residual, 20
                ),
                "recurrence_residual": scientific_text(
                    solved.recurrence_residual, 20
                ),
                "annihilator_coefficients": coefficient_strings,
            }

            next_exponent = start + order + 1
            if next_exponent in exponent_set:
                history = [dyadic_roots[exponent] for exponent in required[1:]]
                predicted = predict_next(
                    history,
                    center=solved.center,
                    coefficients=solved.coefficients,
                )
                observed = dyadic_roots[next_exponent]
                error = predicted - observed
                observed_step = observed - values[-1]
                row["rolling_holdout"] = {
                    "cutoff": 2**next_exponent,
                    "prediction": mp_text(predicted, display_digits),
                    "observed": mp_text(observed, display_digits),
                    "signed_error": scientific_text(error, 20),
                    "absolute_error": scientific_text(abs(error), 20),
                    "relative_to_observed_step": (
                        scientific_text(abs(error / observed_step), 16)
                        if observed_step != 0
                        else None
                    ),
                    "float64_ulp_at_observed": repr(math.ulp(float(observed))),
                    "error_in_float64_ulps": scientific_text(
                        abs(error) / math.ulp(float(observed)), 16
                    ),
                }
            else:
                row["rolling_holdout"] = None
            windows.append(row)

    report: dict[str, Any] = {
        "linear_layers": linear_layers,
        "tail_exponents": [mp_text(mp.mpf(decay) + index, 20) for index in range(linear_layers)],
        "recurrence_order": order,
        "required_consecutive_cutoffs": order + 1,
        "windows": windows,
        "consensus": None,
        "holdout_summary": None,
    }
    if not windows:
        return report

    selected = windows[-max(1, consensus_windows) :]
    centers = [mp.mpf(row["center"]) for row in selected]
    consensus = median_mpf(centers)
    spread = max(centers) - min(centers)
    holdouts = [
        row["rolling_holdout"]
        for row in windows
        if row["rolling_holdout"] is not None
    ]
    holdout_errors = [mp.mpf(item["absolute_error"]) for item in holdouts]
    recent_errors = holdout_errors[-max(1, recent_holdouts) :]
    recent_max = max(recent_errors) if recent_errors else mp.mpf("0")
    empirical_radius = max(spread, recent_max)
    places = safe_decimal_places(empirical_radius, cap=dps - 10)
    integer_digits = max(1, int(mp.floor(mp.log10(abs(consensus)))) + 1) if consensus else 1

    report["consensus"] = {
        "method": "median of latest phase-consistent windows",
        "window_count": len(selected),
        "window_cutoffs": [row["cutoffs"] for row in selected],
        "center": mp_text(consensus, display_digits),
        "center_spread": scientific_text(spread, 20),
        "recent_holdout_count": len(recent_errors),
        "recent_holdout_max_abs_error": scientific_text(recent_max, 20),
        "empirical_radius": scientific_text(empirical_radius, 20),
        "empirical_radius_is_rigorous": False,
        "safe_decimal_places_from_empirical_radius": places,
        "safe_significant_digits_from_empirical_radius": (
            integer_digits + places if places is not None else None
        ),
        "recommended_rounded_center": (
            fixed_decimal(consensus, places) if places is not None else None
        ),
    }

    if holdouts:
        absolute_errors = [mp.mpf(item["absolute_error"]) for item in holdouts]
        report["holdout_summary"] = {
            "count": len(holdouts),
            "max_abs_error": scientific_text(max(absolute_errors), 20),
            "median_abs_error": scientific_text(median_mpf(absolute_errors), 20),
            "rms_abs_error": scientific_text(
                mp.sqrt(mp.fsum(error * error for error in absolute_errors) / len(absolute_errors)),
                20,
            ),
            "latest_abs_error": scientific_text(absolute_errors[-1], 20),
        }
    return report


def generate_predictions(
    dyadic_roots: dict[int, Any],
    *,
    center: Any,
    decay: Any,
    linear_layers: int,
    count: int,
    dps: int,
    display_digits: int,
) -> list[dict[str, Any]]:
    if count <= 0:
        return []
    order = 2 * linear_layers
    exponents = sorted(dyadic_roots)
    if len(exponents) < order:
        return []
    last = exponents[-1]
    required = list(range(last - order + 1, last + 1))
    if not all(exponent in dyadic_roots for exponent in required):
        return []
    history = [dyadic_roots[exponent] for exponent in required]
    coefficients = annihilator_coefficients(center, decay, linear_layers)
    rows: list[dict[str, Any]] = []
    current_exponent = last
    for _ in range(count):
        prediction = predict_next(
            history,
            center=center,
            coefficients=coefficients,
        )
        current_exponent += 1
        rows.append(
            {
                "cutoff": 2**current_exponent,
                "exponent": current_exponent,
                "predicted_root": mp_text(prediction, display_digits),
                "predicted_drift_from_center": scientific_text(
                    prediction - center, 20
                ),
                "prediction_is_certified": False,
            }
        )
        history = history[1:] + [prediction]
    return rows


def layer_label(layer_count: int) -> str:
    if layer_count == 1:
        return "leading linear tail"
    if layer_count == 2:
        return "canonical two-layer linear tail"
    return f"exploratory {layer_count}-layer linear tail"


def markdown_report(payload: dict[str, Any]) -> str:
    canonical = payload["canonical_result"]
    lines: list[str] = []
    lines.append("# Native-Carry Asymptotic Tail Corrector")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(f"- Camera: `{payload['operator_context'].get('camera')}`")
    lines.append(f"- Geometry: `{payload['operator_context'].get('geometry')}`")
    lines.append(f"- Dyadic cutoff points used: `{payload['data_audit']['dyadic_point_count']}`")
    lines.append(f"- Canonical model: `{canonical['label']}`")
    consensus = canonical.get("consensus")
    if consensus:
        lines.append(
            f"- Corrected center: **`{consensus['center']}`**"
        )
        lines.append(
            f"- Recommended empirical rounding: **`{consensus['recommended_rounded_center']}`**"
        )
        lines.append(
            f"- Empirical model radius: `{consensus['empirical_radius']}` (not a rigorous interval)"
        )
    lines.append("")
    lines.append("The finite operator was not changed. The corrector acts only on the sequence of finite-cutoff minima.")
    lines.append("")
    lines.append("## Native derivation")
    lines.append("")
    lines.append("The native real rotation and centered second difference produce the tail layers")
    lines.append("")
    lines.append("```text")
    lines.append("M^(-(3/2+j)) Rot(-t log M) v_j(t),  j=0,1,2,...")
    lines.append("```")
    lines.append("")
    lines.append("Doubling the cutoff contracts layer `j` by `2^(-(3/2+j))` and rotates it by `-t log 2`. The real sequence is therefore annihilated by a second-order polynomial per layer. The canonical corrector multiplies the first two polynomials and solves the resulting five-point relation phase-consistently.")
    lines.append("")
    lines.append("## Rolling holdouts")
    lines.append("")
    lines.append("| Training cutoffs | Predicted cutoff | Prediction | Observed | Absolute error | Error / float64 ULP |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for window in canonical.get("windows", []):
        holdout = window.get("rolling_holdout")
        if not holdout:
            continue
        cutoffs = ", ".join(str(value) for value in window["cutoffs"])
        lines.append(
            "| "
            + cutoffs
            + f" | {holdout['cutoff']} | `{holdout['prediction']}` | `{holdout['observed']}` | `{holdout['absolute_error']}` | `{holdout['error_in_float64_ulps']}` |"
        )
    if not any(window.get("rolling_holdout") for window in canonical.get("windows", [])):
        lines.append("| — | — | — | — | no holdout available | — |")
    lines.append("")
    lines.append("## Improvement over the leading-layer corrector")
    lines.append("")
    lines.append("| Cutoff | One-layer error | Canonical error | Improvement factor |")
    lines.append("|---:|---:|---:|---:|")
    comparisons = payload.get("model_comparison", {}).get("matched_rolling_holdouts", [])
    for row in comparisons:
        lines.append(
            f"| {row['cutoff']} | `{row['leading_layer_abs_error']}` | `{row['canonical_abs_error']}` | `{row['improvement_factor']}` |"
        )
    if not comparisons:
        lines.append("| — | — | — | no matched holdout |")
    lines.append("")
    lines.append("## Local corrected centers")
    lines.append("")
    lines.append("| Linear layers | Cutoffs | Corrected center | Recurrence residual |")
    lines.append("|---:|---|---:|---:|")
    for report in payload["layer_reports"]:
        for window in report["windows"]:
            cutoffs = ", ".join(str(value) for value in window["cutoffs"])
            lines.append(
                f"| {report['linear_layers']} | {cutoffs} | `{window['center']}` | `{window['recurrence_residual']}` |"
            )
    lines.append("")
    lines.append("## Next model predictions")
    lines.append("")
    lines.append("| Cutoff | Predicted finite minimum | Drift from corrected center |")
    lines.append("|---:|---:|---:|")
    for row in payload.get("predictions", []):
        lines.append(
            f"| {row['cutoff']} | `{row['predicted_root']}` | `{row['predicted_drift_from_center']}` |"
        )
    if not payload.get("predictions"):
        lines.append("| — | no prediction | — |")
    lines.append("")
    lines.append("## Logical status")
    lines.append("")
    lines.append("- Finite-input witnesses: high-precision stationary minima emitted by the unchanged native-carry operator.")
    lines.append("- Derived asymptotic structure: centered second-difference tail and dyadic contraction/rotation.")
    lines.append("- Empirical validation: rolling cutoff holdouts.")
    lines.append("- Not claimed: an interval-certified limit or a rigorous bound for the omitted nonlinear remainder.")
    lines.append("")
    lines.append("## Data audit")
    lines.append("")
    lines.append(f"- All cutoff points: `{payload['data_audit']['all_point_count']}`")
    lines.append(f"- Dyadic points: `{payload['data_audit']['dyadic_point_count']}`")
    lines.append(f"- Non-dyadic points excluded from recurrence: `{payload['data_audit']['non_dyadic_cutoffs']}`")
    lines.append(f"- Missing dyadic cutoffs inside observed span: `{payload['data_audit']['missing_dyadic_cutoffs']}`")
    lines.append(f"- Conflicting duplicate witnesses: `{len(payload['data_audit']['conflicts'])}`")
    lines.append("")
    return "\n".join(lines) + "\n"


def run_self_test() -> int:
    print("SELF-TEST NATIVE-CARRY ASYMPTOTIC CORRECTOR")
    with mp.workdps(100):
        true_center = mp.mpf("12.345678901234567890123456789")
        decay = mp.mpf("1.5")
        theta = true_center * mp.log(2)

        def synthetic(exponent: int, layers: int) -> Any:
            total = true_center
            amplitudes = [
                mp.mpc("3.2", "-1.7"),
                mp.mpc("-5.1", "2.4"),
                mp.mpc("0.7", "0.9"),
            ]
            for layer in range(layers):
                rho = mp.power(2, -(decay + layer))
                mode = amplitudes[layer] * mp.power(
                    rho * mp.e ** (-1j * theta), exponent
                )
                total += mp.re(mode)
            return total

        one_values = [synthetic(exponent, 1) for exponent in range(8, 11)]
        one = phase_consistent_center(
            one_values,
            decay=decay,
            linear_layers=1,
            tolerance=mp.mpf("1e-80"),
        )
        if abs(one.center - true_center) > mp.mpf("1e-70"):
            raise AssertionError("one-layer center recovery failed")

        two_values = [synthetic(exponent, 2) for exponent in range(8, 13)]
        two = phase_consistent_center(
            two_values,
            decay=decay,
            linear_layers=2,
            tolerance=mp.mpf("1e-80"),
        )
        if abs(two.center - true_center) > mp.mpf("1e-70"):
            raise AssertionError("two-layer center recovery failed")
        predicted = predict_next(
            two_values[1:],
            center=two.center,
            coefficients=two.coefficients,
        )
        if abs(predicted - synthetic(13, 2)) > mp.mpf("1e-70"):
            raise AssertionError("two-layer prediction failed")

    print("  PASS one-layer phase-consistent center")
    print("  PASS two-layer annihilator and holdout prediction")
    print("  PASS no Newton and no interval arithmetic")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Correct finite native-carry minima by annihilating damped dyadic "
            "tail rotations; reads precision-ladder JSON ledgers"
        ),
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        type=Path,
        help="precision-ladder JSON files; duplicate cutoffs are merged",
    )
    parser.add_argument(
        "--point",
        action="append",
        default=[],
        help="manual CUTOFF:ROOT witness; may be repeated",
    )
    parser.add_argument("--camera", type=int)
    parser.add_argument(
        "--decay",
        default=DEFAULT_DECAY,
        help="leading real decay exponent; 3/2 follows from the native amplitude and second differences",
    )
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument(
        "--display-digits",
        type=int,
        default=50,
        help="significant digits printed for derived model values; empirical rounding is reported separately",
    )
    parser.add_argument(
        "--max-linear-layers",
        type=int,
        default=3,
        help="report 1..N linear tail layers; layers >2 are exploratory",
    )
    parser.add_argument(
        "--canonical-layers",
        type=int,
        default=2,
        help="layer count used for the published corrected center and predictions",
    )
    parser.add_argument("--consensus-windows", type=int, default=3)
    parser.add_argument("--recent-holdouts", type=int, default=2)
    parser.add_argument("--predict-doublings", type=int, default=3)
    parser.add_argument("--json-out", type=Path)
    parser.add_argument("--markdown-out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.self_test:
        return run_self_test()

    if not args.inputs and not args.point:
        raise SystemExit("provide at least one JSON input or --point CUTOFF:ROOT")
    if args.camera is not None and args.camera < 2:
        raise SystemExit("camera must be >= 2")
    if args.dps < 50:
        raise SystemExit("dps must be >= 50")
    if not 17 <= args.display_digits <= args.dps:
        raise SystemExit("display-digits must lie between 17 and dps")
    if args.max_linear_layers < 1:
        raise SystemExit("max-linear-layers must be >= 1")
    if not 1 <= args.canonical_layers <= args.max_linear_layers:
        raise SystemExit("canonical-layers must lie in 1..max-linear-layers")
    if args.consensus_windows < 1 or args.recent_holdouts < 1:
        raise SystemExit("consensus-windows and recent-holdouts must be positive")
    if args.predict_doublings < 0:
        raise SystemExit("predict-doublings must be nonnegative")

    try:
        decay = mp.mpf(args.decay)
        if decay <= 0:
            raise ValueError
    except Exception as exc:
        raise SystemExit("decay must be a positive real number") from exc

    with mp.workdps(args.dps):
        try:
            points, conflicts, source_metadata = load_points(
                args.inputs,
                args.point,
                requested_camera=args.camera,
            )
        except Exception as exc:
            raise SystemExit(str(exc)) from exc

        all_cutoffs = sorted(points)
        dyadic_points = {
            power_of_two_exponent(cutoff): mp.mpf(point.root_text)
            for cutoff, point in points.items()
            if is_power_of_two(cutoff)
        }
        non_dyadic = [cutoff for cutoff in all_cutoffs if not is_power_of_two(cutoff)]
        dyadic_exponents = sorted(dyadic_points)
        missing: list[int] = []
        if dyadic_exponents:
            missing = [
                2**exponent
                for exponent in range(dyadic_exponents[0], dyadic_exponents[-1] + 1)
                if exponent not in dyadic_points
            ]

        layer_reports = [
            build_layer_report(
                dyadic_points,
                decay=decay,
                linear_layers=layers,
                dps=args.dps,
                consensus_windows=args.consensus_windows,
                recent_holdouts=args.recent_holdouts,
                display_digits=args.display_digits,
            )
            for layers in range(1, args.max_linear_layers + 1)
        ]
        canonical = layer_reports[args.canonical_layers - 1]
        if canonical["consensus"] is None:
            required = 2 * args.canonical_layers + 1
            raise SystemExit(
                f"canonical {args.canonical_layers}-layer corrector needs at least "
                f"{required} consecutive dyadic cutoffs"
            )

        canonical_center = mp.mpf(canonical["consensus"]["center"])
        predictions = generate_predictions(
            dyadic_points,
            center=canonical_center,
            decay=decay,
            linear_layers=args.canonical_layers,
            count=args.predict_doublings,
            dps=args.dps,
            display_digits=args.display_digits,
        )

        # Compare the canonical model against the leading-layer model on exactly
        # the same rolling holdouts.  This avoids mixing different target cutoffs.
        holdouts_by_layer: dict[int, dict[int, Any]] = {}
        for report in layer_reports:
            layer_map: dict[int, Any] = {}
            for window in report["windows"]:
                holdout = window.get("rolling_holdout")
                if holdout is not None:
                    layer_map[int(holdout["cutoff"])] = mp.mpf(holdout["absolute_error"])
            holdouts_by_layer[int(report["linear_layers"])] = layer_map
        matched_improvements: list[dict[str, Any]] = []
        leading_map = holdouts_by_layer.get(1, {})
        canonical_map = holdouts_by_layer.get(args.canonical_layers, {})
        for cutoff in sorted(set(leading_map) & set(canonical_map)):
            canonical_error = canonical_map[cutoff]
            leading_error = leading_map[cutoff]
            matched_improvements.append(
                {
                    "cutoff": cutoff,
                    "leading_layer_abs_error": scientific_text(leading_error, 20),
                    "canonical_abs_error": scientific_text(canonical_error, 20),
                    "improvement_factor": (
                        scientific_text(leading_error / canonical_error, 16)
                        if canonical_error != 0
                        else None
                    ),
                }
            )

        point_rows = []
        for cutoff in all_cutoffs:
            point = points[cutoff]
            point_rows.append(
                {
                    "cutoff": cutoff,
                    "dyadic": is_power_of_two(cutoff),
                    "exponent": (
                        power_of_two_exponent(cutoff)
                        if is_power_of_two(cutoff)
                        else None
                    ),
                    "root": point.root_text,
                    "score": point.score_text,
                    "resultant_norm": point.resultant_norm_text,
                    "source": point.source,
                    "root_digits": point.root_digits,
                }
            )

        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "status": "FINITE_NATIVE_CARRY_ASYMPTOTIC_CORRECTION_AUDIT",
            "interpretation": {
                "operator_changed": False,
                "newton": False,
                "interval_arithmetic": False,
                "finite_cutoff_inputs": True,
                "rigorous_remainder_bound": False,
                "uncertainty_kind": "empirical rolling-holdout and window spread",
                "canonical_model": "two linear tail layers" if args.canonical_layers == 2 else layer_label(args.canonical_layers),
                "higher_linear_layers_warning": (
                    "nonlinear root-shift terms may enter before the third linear layer"
                ),
            },
            "operator_context": {
                "camera": source_metadata["detected_camera"],
                "geometry": source_metadata["detected_geometry"],
                "state": source_metadata["detected_state"],
            },
            "native_tail_model": {
                "real_state": "f_t(n)=n^(-1/2)*Rot(-t*log(n))*e1",
                "centered_difference": "f(n-r)-2f(n)+f(n+r)=r^2*f''(n)+r^4*f''''(n)/12+...",
                "real_rotational_tail_layers": "M^(-(3/2+j))*Rot(-t*log(M))*v_j(t), j=0,1,2,...",
                "leading_real_decay": mp_text(decay, 20),
                "dyadic_real_map_layer_j": "2^(-(decay+j))*Rot(-T*log(2))",
                "annihilator_layer_j": "L^2-2*rho_j*cos(T*log(2))*L+rho_j^2",
                "root_drift_assumption": "nondegenerate limiting stationary minimum and dominance of the displayed tail layers",
            },
            "data_audit": {
                "all_point_count": len(points),
                "dyadic_point_count": len(dyadic_points),
                "all_cutoffs": all_cutoffs,
                "dyadic_cutoffs": [2**exponent for exponent in dyadic_exponents],
                "non_dyadic_cutoffs": non_dyadic,
                "missing_dyadic_cutoffs": missing,
                "conflicts": conflicts,
                "points": point_rows,
                "sources": source_metadata["sources"],
            },
            "layer_reports": layer_reports,
            "model_comparison": {
                "baseline": "one leading linear tail layer",
                "canonical": layer_label(args.canonical_layers),
                "matched_rolling_holdouts": matched_improvements,
            },
            "canonical_result": {
                "linear_layers": args.canonical_layers,
                "label": layer_label(args.canonical_layers),
                "consensus": canonical["consensus"],
                "holdout_summary": canonical["holdout_summary"],
                "windows": canonical["windows"],
            },
            "predictions": predictions,
            "runtime": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "mpmath": mp.__version__,
                "dps": args.dps,
                "script_sha256": file_sha256(Path(__file__).resolve()),
            },
        }

        consensus = canonical["consensus"]
        print("=" * 120)
        print(" NATIVE-CARRY ASYMPTOTIC TAIL CORRECTOR — DYADIC ROTATION/CONTRACTION, NO NEWTON")
        print("=" * 120)
        print(f" camera                  : {source_metadata['detected_camera']}")
        print(f" geometry                : {source_metadata['detected_geometry']}")
        print(f" cutoff witnesses        : {len(points)} total, {len(dyadic_points)} dyadic")
        print(f" dyadic span             : {min(2**e for e in dyadic_exponents):,} .. {max(2**e for e in dyadic_exponents):,}")
        print(f" leading decay           : {mp_text(decay, 20)}")
        print(f" canonical linear layers : {args.canonical_layers}")
        print(f" recurrence order        : {canonical['recurrence_order']}")
        print("-" * 120)
        print(" Rolling holdouts:")
        for window in canonical["windows"]:
            holdout = window.get("rolling_holdout")
            if holdout:
                print(
                    f"  train={window['cutoffs']} -> M={holdout['cutoff']:<8d} "
                    f"pred={holdout['prediction']}  error={holdout['signed_error']} "
                    f"({holdout['error_in_float64_ulps']} float64 ULP)"
                )
        print("-" * 120)
        print(f" corrected center        : {consensus['center']}")
        print(f" recent-window spread    : {consensus['center_spread']}")
        print(f" recent holdout max err  : {consensus['recent_holdout_max_abs_error']}")
        print(f" empirical radius        : {consensus['empirical_radius']}  (not an interval proof)")
        print(f" recommended rounding    : {consensus['recommended_rounded_center']}")
        if predictions:
            print("-" * 120)
            print(" Next model predictions:")
            for row in predictions:
                print(
                    f"  M={row['cutoff']:<10d} t={row['predicted_root']} "
                    f"drift={row['predicted_drift_from_center']}"
                )
        if non_dyadic or missing:
            print("-" * 120)
            print(f" non-dyadic excluded     : {non_dyadic or '-'}")
            print(f" missing dyadic cutoffs  : {missing or '-'}")
        print("-" * 120)
        print(" STATUS")
        print(" * finite minima are unchanged operator outputs")
        print(" * canonical correction annihilates the real rotational layers with decays M^(-3/2) and M^(-5/2)")
        print(" * uncertainty is empirical; no rigorous remainder bound is claimed")

        if args.json_out is not None:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            temporary = args.json_out.with_suffix(args.json_out.suffix + ".tmp")
            temporary.write_text(
                json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            temporary.replace(args.json_out)
            print(f" JSON written             : {args.json_out}")
        if args.markdown_out is not None:
            args.markdown_out.parent.mkdir(parents=True, exist_ok=True)
            temporary = args.markdown_out.with_suffix(args.markdown_out.suffix + ".tmp")
            temporary.write_text(markdown_report(payload), encoding="utf-8")
            temporary.replace(args.markdown_out)
            print(f" Markdown written         : {args.markdown_out}")
        print("=" * 120)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
