#!/usr/bin/env python3
"""High-precision laboratory for the native-carry primitive operator over R^2.

Purpose
-------
The existing scanner is deliberately excellent at the *wide* search:
CUDA/CPU float64 evaluates a large decimal grid and identifies deep local
valleys.  This companion script performs the *narrow* inspection after a
candidate is known.

It preserves the same finite primitive operator:

    psi_t(n) = n^(-1/2) (cos(-t log n), sin(-t log n)).

Camera 2:
    seed psi_t(1), centers 4,8,...,4M, radius 1.

Natural camera b >= 3:
    h=floor(b/2), seeds 1,...,h, centers b,2b,...,Mb,
    radii 1,...,h, including the antipodal radius b/2 when b is even.

For each bracket z=(psi(c-r)-2 psi(c)+psi(c+r)),

    R(t) = sum z_e,
    E(t) = sum ||z_e||^2,
    score(t) = ||R(t)||^2 / (N E(t)).

What changes is only the arithmetic used to inspect a small neighborhood:
mpmath arbitrary precision replaces float64.  There is no calibration, no
post-bracket map, no interval arithmetic, and no Newton iteration.

The local minimum is found by a bracketed Ridder solve.  For the canonical
score objective, the solved scalar condition is

    H(t) = F'(t) E(t) - F(t) E'(t) = 0,
    F(t) = ||R(t)||^2,

whose sign is the sign of score'(t), since N*E(t)^2 > 0.  Ridder's method is
bracketed and derivative-free as a root solver.  R'(t) and E'(t) are evaluated
analytically from the unchanged rotation merely to form H(t).

Important interpretation
------------------------
This script locates a minimum of a *finite-cutoff* operator.  Arbitrary
arithmetic can determine that finite minimum to many digits, but the digits
that survive a cutoff ladder are the meaningful convergence diagnostic.
A tiny score is not itself a count of correct digits in t.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import platform
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

try:
    import mpmath as mp
except Exception as exc:  # pragma: no cover - dependency error path
    raise SystemExit(
        "mpmath is required. Install it with: python3 -m pip install mpmath"
    ) from exc


SCHEMA = "org.native-carry.primitive-r2-precision-ladder/v1"
DEFAULT_DPS = (40, 70, 110)
DEFAULT_OPERATOR_NAME = "native_carry_primitive_real_operator_all_bases.py"


@dataclass(frozen=True)
class CameraGeometry:
    camera: int
    cutoff: int
    half_range: int
    geometry: str
    seeds: tuple[int, ...]
    brackets: tuple[tuple[int, int, int], ...]
    unique_n: tuple[int, ...]
    seed_indices: tuple[int, ...]
    bracket_indices: tuple[tuple[int, int, int], ...]
    coefficients: tuple[int, ...]
    coordinate_count: int
    largest_center: int


@dataclass(frozen=True)
class EvalRecord:
    t: Any
    resultant_x: Any
    resultant_y: Any
    derivative_x: Any
    derivative_y: Any
    second_x: Any | None
    second_y: Any | None
    resultant_norm_sq: Any
    resultant_derivative: Any
    energy: Any | None
    energy_derivative: Any | None
    score: Any | None
    score_stationary_numerator: Any | None


@dataclass(frozen=True)
class RootResult:
    root: Any
    left: Any
    right: Any
    f_root: Any
    iterations: int
    evaluations: int
    estimated_error: Any


@dataclass(frozen=True)
class CoarseCandidate:
    t_decimal: str
    score: float
    grid_index: int
    left_score: float | None
    right_score: float | None
    backend: str
    elapsed_seconds: float
    grid_step: str


class PreparedCamera:
    """Camera constants prepared at the active mpmath precision."""

    def __init__(self, geometry: CameraGeometry) -> None:
        self.geometry = geometry
        self.logs = [mp.log(value) for value in geometry.unique_n]
        self.amplitudes = [1 / mp.sqrt(value) for value in geometry.unique_n]
        self.sparse_amplitudes = [
            mp.mpf(coefficient) * amplitude
            for coefficient, amplitude in zip(
                geometry.coefficients, self.amplitudes, strict=True
            )
        ]

    def evaluate_resultant(
        self,
        t: Any,
        *,
        second: bool = False,
    ) -> EvalRecord:
        """Evaluate R, R', and optionally R'' by exact coefficient regrouping."""
        t = mp.mpf(t)
        rx = mp.mpf("0")
        ry = mp.mpf("0")
        dx = mp.mpf("0")
        dy = mp.mpf("0")
        ddx = mp.mpf("0") if second else None
        ddy = mp.mpf("0") if second else None

        for amplitude, log_value in zip(
            self.sparse_amplitudes, self.logs, strict=True
        ):
            if not amplitude:
                continue
            cosine, sine = mp.cos_sin(-t * log_value)
            x = amplitude * cosine
            y = amplitude * sine
            rx += x
            ry += y
            dx += amplitude * log_value * sine
            dy -= amplitude * log_value * cosine
            if second:
                log_sq = log_value * log_value
                assert ddx is not None and ddy is not None
                ddx -= amplitude * log_sq * cosine
                ddy -= amplitude * log_sq * sine

        norm_sq = rx * rx + ry * ry
        f_prime = 2 * (rx * dx + ry * dy)
        return EvalRecord(
            t=t,
            resultant_x=rx,
            resultant_y=ry,
            derivative_x=dx,
            derivative_y=dy,
            second_x=ddx,
            second_y=ddy,
            resultant_norm_sq=norm_sq,
            resultant_derivative=f_prime,
            energy=None,
            energy_derivative=None,
            score=None,
            score_stationary_numerator=None,
        )

    def evaluate_score(
        self,
        t: Any,
        *,
        second: bool = False,
    ) -> EvalRecord:
        """Evaluate the unchanged finite operator and score derivative in R^2.

        Every unique psi_t(n) is rotated once.  Brackets and their energies are
        then assembled from index triples.  The resultant itself is accumulated
        with the exact integer coefficient map, which is only a reassociation of
        the same finite sum.
        """
        t = mp.mpf(t)
        count = len(self.logs)
        x_values: list[Any] = [None] * count
        y_values: list[Any] = [None] * count
        dx_values: list[Any] = [None] * count
        dy_values: list[Any] = [None] * count

        rx = mp.mpf("0")
        ry = mp.mpf("0")
        dx = mp.mpf("0")
        dy = mp.mpf("0")
        ddx = mp.mpf("0") if second else None
        ddy = mp.mpf("0") if second else None

        for index, (amplitude, log_value, coefficient) in enumerate(
            zip(
                self.amplitudes,
                self.logs,
                self.geometry.coefficients,
                strict=True,
            )
        ):
            cosine, sine = mp.cos_sin(-t * log_value)
            x = amplitude * cosine
            y = amplitude * sine
            x_prime = amplitude * log_value * sine
            y_prime = -amplitude * log_value * cosine
            x_values[index] = x
            y_values[index] = y
            dx_values[index] = x_prime
            dy_values[index] = y_prime

            if coefficient:
                rx += coefficient * x
                ry += coefficient * y
                dx += coefficient * x_prime
                dy += coefficient * y_prime
                if second:
                    assert ddx is not None and ddy is not None
                    log_sq = log_value * log_value
                    ddx -= coefficient * log_sq * x
                    ddy -= coefficient * log_sq * y

        # Seed rotations have exactly constant norm 1/n, hence zero E'.
        energy = mp.fsum(
            self.amplitudes[index] * self.amplitudes[index]
            for index in self.geometry.seed_indices
        )
        energy_prime = mp.mpf("0")

        for left, center, right in self.geometry.bracket_indices:
            bracket_x = x_values[left] - 2 * x_values[center] + x_values[right]
            bracket_y = y_values[left] - 2 * y_values[center] + y_values[right]
            bracket_dx = (
                dx_values[left] - 2 * dx_values[center] + dx_values[right]
            )
            bracket_dy = (
                dy_values[left] - 2 * dy_values[center] + dy_values[right]
            )
            energy += bracket_x * bracket_x + bracket_y * bracket_y
            energy_prime += 2 * (
                bracket_x * bracket_dx + bracket_y * bracket_dy
            )

        norm_sq = rx * rx + ry * ry
        f_prime = 2 * (rx * dx + ry * dy)
        score = norm_sq / (mp.mpf(self.geometry.coordinate_count) * energy)
        stationary = f_prime * energy - norm_sq * energy_prime
        return EvalRecord(
            t=t,
            resultant_x=rx,
            resultant_y=ry,
            derivative_x=dx,
            derivative_y=dy,
            second_x=ddx,
            second_y=ddy,
            resultant_norm_sq=norm_sq,
            resultant_derivative=f_prime,
            energy=energy,
            energy_derivative=energy_prime,
            score=score,
            score_stationary_numerator=stationary,
        )


def build_geometry(camera: int, cutoff: int) -> CameraGeometry:
    camera = int(camera)
    cutoff = int(cutoff)
    if camera < 2:
        raise ValueError("camera must be >= 2")
    if cutoff < 1:
        raise ValueError("cutoff must be >= 1")

    if camera == 2:
        half_range = 1
        geometry_name = "c2_aligned_centers_4m"
        seeds = (1,)
        brackets = tuple((4 * m - 1, 4 * m, 4 * m + 1) for m in range(1, cutoff + 1))
        largest_center = 4 * cutoff
    else:
        half_range = camera // 2
        geometry_name = (
            "natural_saturated_even_antipode"
            if camera % 2 == 0
            else "natural_saturated_odd_width"
        )
        seeds = tuple(range(1, half_range + 1))
        brackets = tuple(
            (camera * m - radius, camera * m, camera * m + radius)
            for m in range(1, cutoff + 1)
            for radius in range(1, half_range + 1)
        )
        largest_center = camera * cutoff

    coefficient_map: dict[int, int] = {}

    def add_coefficient(n: int, coefficient: int) -> None:
        coefficient_map[n] = coefficient_map.get(n, 0) + coefficient

    for n in seeds:
        add_coefficient(n, 1)
    for left, center, right in brackets:
        add_coefficient(left, 1)
        add_coefficient(center, -2)
        add_coefficient(right, 1)

    unique_n = tuple(
        sorted(
            set(seeds).union(
                n for bracket in brackets for n in bracket
            )
        )
    )
    index_by_n = {n: index for index, n in enumerate(unique_n)}
    seed_indices = tuple(index_by_n[n] for n in seeds)
    bracket_indices = tuple(
        (index_by_n[left], index_by_n[center], index_by_n[right])
        for left, center, right in brackets
    )
    coefficients = tuple(coefficient_map.get(n, 0) for n in unique_n)

    return CameraGeometry(
        camera=camera,
        cutoff=cutoff,
        half_range=half_range,
        geometry=geometry_name,
        seeds=seeds,
        brackets=brackets,
        unique_n=unique_n,
        seed_indices=seed_indices,
        bracket_indices=bracket_indices,
        coefficients=coefficients,
        coordinate_count=len(seeds) + len(brackets),
        largest_center=largest_center,
    )


def direct_resultant(prepared: PreparedCamera, t: Any) -> tuple[Any, Any]:
    """Literal seed+bracket sum, used only by self-tests."""
    t = mp.mpf(t)

    def state(index: int) -> tuple[Any, Any]:
        amplitude = prepared.amplitudes[index]
        cosine, sine = mp.cos_sin(-t * prepared.logs[index])
        return amplitude * cosine, amplitude * sine

    rx = mp.mpf("0")
    ry = mp.mpf("0")
    for index in prepared.geometry.seed_indices:
        x, y = state(index)
        rx += x
        ry += y
    for left, center, right in prepared.geometry.bracket_indices:
        lx, ly = state(left)
        cx, cy = state(center)
        rx_, ry_ = state(right)
        rx += lx - 2 * cx + rx_
        ry += ly - 2 * cy + ry_
    return rx, ry


class ObjectiveEvaluator:
    def __init__(self, prepared: PreparedCamera, objective: str, key_digits: int) -> None:
        self.prepared = prepared
        self.objective = objective
        self.key_digits = key_digits
        self.cache: dict[str, EvalRecord] = {}

    def _key(self, t: Any) -> str:
        return mp.nstr(mp.mpf(t), n=self.key_digits, strip_zeros=False)

    def evaluate(self, t: Any, *, second: bool = False) -> EvalRecord:
        # Second derivatives are requested only for the final diagnostic and are
        # intentionally not mixed into the root-evaluation cache.
        if second:
            if self.objective == "score":
                return self.prepared.evaluate_score(t, second=True)
            return self.prepared.evaluate_resultant(t, second=True)

        key = self._key(t)
        record = self.cache.get(key)
        if record is None:
            if self.objective == "score":
                record = self.prepared.evaluate_score(t, second=False)
            else:
                record = self.prepared.evaluate_resultant(t, second=False)
            self.cache[key] = record
        return record

    def stationary(self, t: Any) -> Any:
        record = self.evaluate(t)
        if self.objective == "score":
            assert record.score_stationary_numerator is not None
            return record.score_stationary_numerator
        # Half of d||R||^2/dt; the factor 2 does not change the root.
        return record.resultant_derivative / 2

    def objective_value(self, t: Any) -> Any:
        record = self.evaluate(t)
        if self.objective == "score":
            assert record.score is not None
            return record.score
        return record.resultant_norm_sq


def _negative_to_positive(left_value: Any, right_value: Any) -> bool:
    return left_value <= 0 and right_value >= 0 and not (
        left_value == 0 and right_value == 0
    )


def locate_minimum_bracket(
    evaluator: ObjectiveEvaluator,
    center: Any,
    half_width: Any,
    *,
    max_expansions: int,
    fallback_samples: int = 16,
) -> tuple[Any, Any, Any, Any, int]:
    """Find a negative-to-positive stationary crossing near center."""
    center = mp.mpf(center)
    width = abs(mp.mpf(half_width))
    if width == 0:
        raise ValueError("half-width must be positive")

    for expansion in range(max_expansions + 1):
        left = center - width
        right = center + width
        f_left = evaluator.stationary(left)
        f_right = evaluator.stationary(right)
        if _negative_to_positive(f_left, f_right):
            return left, right, f_left, f_right, expansion

        # Fallback sampling is used only when the two endpoints do not bracket.
        # It selects the closest negative-to-positive crossing, excluding maxima.
        step = (right - left) / fallback_samples
        points = [left + index * step for index in range(fallback_samples + 1)]
        values = [evaluator.stationary(point) for point in points]
        candidates: list[tuple[Any, int]] = []
        for index in range(fallback_samples):
            if _negative_to_positive(values[index], values[index + 1]):
                midpoint = (points[index] + points[index + 1]) / 2
                candidates.append((abs(midpoint - center), index))
        if candidates:
            _, index = min(candidates, key=lambda item: item[0])
            return (
                points[index],
                points[index + 1],
                values[index],
                values[index + 1],
                expansion,
            )
        width *= 2

    raise RuntimeError(
        "could not bracket a local minimum; enlarge --half-width or provide a better coarse candidate"
    )


def ridder_root(
    function: Callable[[Any], Any],
    left: Any,
    right: Any,
    f_left: Any,
    f_right: Any,
    *,
    x_tolerance: Any,
    max_steps: int,
    evaluation_count: Callable[[], int],
) -> RootResult:
    """Bracketed Ridder root solve, with bisection fallbacks and no Newton step."""
    left = mp.mpf(left)
    right = mp.mpf(right)
    f_left = mp.mpf(f_left)
    f_right = mp.mpf(f_right)
    x_tolerance = abs(mp.mpf(x_tolerance))

    if left > right:
        left, right = right, left
        f_left, f_right = f_right, f_left
    if f_left == 0:
        return RootResult(left, left, left, f_left, 0, evaluation_count(), mp.mpf("0"))
    if f_right == 0:
        return RootResult(right, right, right, f_right, 0, evaluation_count(), mp.mpf("0"))
    if f_left * f_right > 0:
        raise ValueError("Ridder requires a sign-changing bracket")

    best_x = left if abs(f_left) <= abs(f_right) else right
    best_f = f_left if abs(f_left) <= abs(f_right) else f_right
    previous_candidate: Any | None = None

    for iteration in range(1, max_steps + 1):
        midpoint = (left + right) / 2
        f_mid = function(midpoint)
        if abs(f_mid) < abs(best_f):
            best_x, best_f = midpoint, f_mid

        radicand = f_mid * f_mid - f_left * f_right
        if radicand <= 0:
            candidate = midpoint
        else:
            denominator = mp.sqrt(radicand)
            sign = mp.mpf("1") if f_left >= f_right else mp.mpf("-1")
            candidate = midpoint + (midpoint - left) * sign * f_mid / denominator
            if not (left < candidate < right):
                candidate = midpoint

        f_candidate = function(candidate)
        if abs(f_candidate) < abs(best_f):
            best_x, best_f = candidate, f_candidate

        step_error = (
            abs(candidate - previous_candidate)
            if previous_candidate is not None
            else abs(right - left)
        )
        if f_candidate == 0 or (
            previous_candidate is not None and step_error <= x_tolerance
        ):
            return RootResult(
                candidate,
                left,
                right,
                f_candidate,
                iteration,
                evaluation_count(),
                step_error,
            )
        previous_candidate = candidate

        if f_mid * f_candidate < 0:
            left, f_left = midpoint, f_mid
            right, f_right = candidate, f_candidate
        elif f_left * f_candidate < 0:
            right, f_right = candidate, f_candidate
        elif f_right * f_candidate < 0:
            left, f_left = candidate, f_candidate
        else:
            # Degenerate arithmetic case: keep the half containing the sign change.
            if f_left * f_mid <= 0:
                right, f_right = midpoint, f_mid
            else:
                left, f_left = midpoint, f_mid

        if left > right:
            left, right = right, left
            f_left, f_right = f_right, f_left

        if abs(right - left) <= x_tolerance:
            return RootResult(
                best_x,
                left,
                right,
                best_f,
                iteration,
                evaluation_count(),
                abs(right - left),
            )

    raise RuntimeError(
        f"Ridder did not converge in {max_steps} steps; final width={mp.nstr(abs(right-left), 12)}"
    )


def parse_int_list(text: str, *, name: str, minimum: int = 1) -> tuple[int, ...]:
    try:
        values = tuple(int(piece.strip()) for piece in text.split(",") if piece.strip())
    except ValueError as exc:
        raise ValueError(f"invalid {name}: {exc}") from exc
    if not values:
        raise ValueError(f"{name} cannot be empty")
    if any(value < minimum for value in values):
        raise ValueError(f"every {name} value must be >= {minimum}")
    return values


def significant_agreement(first: Any, second: Any, cap: int) -> int:
    first = mp.mpf(first)
    second = mp.mpf(second)
    difference = abs(first - second)
    if difference == 0:
        return cap
    scale = max(mp.mpf("1"), abs(first), abs(second))
    digits = int(mp.floor(-mp.log10(difference / scale)))
    return max(0, min(cap, digits))


def decimal_place_agreement(first: Any, second: Any, cap: int) -> int:
    difference = abs(mp.mpf(first) - mp.mpf(second))
    if difference == 0:
        return cap
    digits = int(mp.floor(-mp.log10(difference)))
    return max(0, min(cap, digits))


def mp_text(value: Any, digits: int) -> str:
    return mp.nstr(value, n=max(8, digits), strip_zeros=False)


def scientific_text(value: Any, digits: int = 18) -> str:
    return mp.nstr(value, n=max(8, digits), min_fixed=0, max_fixed=0)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def import_operator(path: Path) -> Any:
    if not path.exists():
        raise FileNotFoundError(path)
    spec = importlib.util.spec_from_file_location("native_carry_float64_operator", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not import operator from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def discover_coarse_candidate(args: argparse.Namespace) -> CoarseCandidate:
    operator_path = Path(args.operator).expanduser().resolve()
    module = import_operator(operator_path)
    grid = module.DecimalGrid.from_strings(
        args.coarse_tmin, args.coarse_tmax, args.coarse_grid
    )
    times = grid.float_values()
    model = module.build_camera_model(args.camera, args.cutoff)
    started = time.perf_counter()
    scores, backend_name, _ = module.scan_camera(
        times,
        model,
        args.coarse_backend,
        args.workers,
        args.cpu_batch,
        args.state_block,
        args.gpu_batch,
        args.gpu_threads,
    )
    indices, prominence_by_index = module.detect_grid_candidates(
        scores, args.prominence, args.prominence_window
    )
    candidates = module.canonicalize_candidates_on_cpu(
        indices,
        times,
        model,
        args.state_block,
        prominence_by_index,
    )
    if not candidates:
        raise RuntimeError("the coarse window contains no local grid minimum")
    best = min(candidates, key=lambda row: row[1])
    index, score, left_score, right_score, _ = best
    return CoarseCandidate(
        t_decimal=grid.text_at(int(index)),
        score=float(score),
        grid_index=int(index),
        left_score=left_score,
        right_score=right_score,
        backend=backend_name,
        elapsed_seconds=time.perf_counter() - started,
        grid_step=str(args.coarse_grid),
    )


def stage_payload(
    *,
    geometry: CameraGeometry,
    requested_dps: int,
    work_dps: int,
    objective: str,
    center: Any,
    initial_half_width: Any,
    root_result: RootResult,
    bracket_expansions: int,
    final_record: EvalRecord,
    elapsed_seconds: float,
    stable_significant_digits: int | None,
    stable_decimal_places: int | None,
) -> dict[str, Any]:
    root = root_result.root
    residual_norm = mp.sqrt(final_record.resultant_norm_sq)
    derivative_norm = mp.sqrt(
        final_record.derivative_x * final_record.derivative_x
        + final_record.derivative_y * final_record.derivative_y
    )
    raw_curvature = None
    if final_record.second_x is not None and final_record.second_y is not None:
        raw_curvature = 2 * (
            derivative_norm * derivative_norm
            + final_record.resultant_x * final_record.second_x
            + final_record.resultant_y * final_record.second_y
        )

    float_root = float(root)
    float_ulp = math.ulp(float_root)
    root_minus_float = root - mp.mpf(float_root)
    score_digits = None
    if final_record.score is not None and final_record.score > 0:
        score_digits = float(-mp.log10(final_record.score))

    return {
        "camera": geometry.camera,
        "cutoff": geometry.cutoff,
        "geometry": geometry.geometry,
        "coordinate_count": geometry.coordinate_count,
        "unique_state_count": len(geometry.unique_n),
        "largest_center": geometry.largest_center,
        "objective": objective,
        "requested_dps": requested_dps,
        "work_dps": work_dps,
        "search_center": mp_text(center, requested_dps + 8),
        "initial_half_width": mp_text(initial_half_width, requested_dps + 8),
        "root": mp_text(root, requested_dps + 8),
        "root_bracket_left": mp_text(root_result.left, requested_dps + 8),
        "root_bracket_right": mp_text(root_result.right, requested_dps + 8),
        "root_bracket_width": scientific_text(
            abs(root_result.right - root_result.left), 18
        ),
        "root_error_estimate": scientific_text(root_result.estimated_error, 18),
        "stationary_numerator": scientific_text(root_result.f_root, 18),
        "ridder_iterations": root_result.iterations,
        "function_evaluations": root_result.evaluations,
        "bracket_expansions": bracket_expansions,
        "resultant": [
            scientific_text(final_record.resultant_x, requested_dps),
            scientific_text(final_record.resultant_y, requested_dps),
        ],
        "resultant_norm": scientific_text(residual_norm, requested_dps),
        "resultant_norm_sq": scientific_text(
            final_record.resultant_norm_sq, requested_dps
        ),
        "resultant_derivative_norm": scientific_text(
            derivative_norm, requested_dps
        ),
        "raw_resultant_curvature": (
            scientific_text(raw_curvature, requested_dps)
            if raw_curvature is not None
            else None
        ),
        "energy": (
            scientific_text(final_record.energy, requested_dps)
            if final_record.energy is not None
            else None
        ),
        "energy_derivative": (
            scientific_text(final_record.energy_derivative, requested_dps)
            if final_record.energy_derivative is not None
            else None
        ),
        "score": (
            scientific_text(final_record.score, requested_dps)
            if final_record.score is not None
            else None
        ),
        "minus_log10_score": score_digits,
        "float64_projection": {
            "root": repr(float_root),
            "ulp_at_root": repr(float_ulp),
            "high_precision_minus_float64": scientific_text(root_minus_float, 18),
        },
        "agreement_with_previous_stage": {
            "significant_digits": stable_significant_digits,
            "decimal_places": stable_decimal_places,
        },
        "elapsed_seconds": elapsed_seconds,
    }


def run_precision_stage(
    *,
    geometry: CameraGeometry,
    center_text: str,
    half_width_text: str,
    requested_dps: int,
    guard_digits: int,
    objective: str,
    max_steps: int,
    max_expansions: int,
    previous_root_text: str | None,
    previous_dps: int | None,
    previous_width_text: str | None,
) -> tuple[dict[str, Any], str, str]:
    work_dps = requested_dps + guard_digits
    started = time.perf_counter()
    with mp.workdps(work_dps):
        prepared = PreparedCamera(geometry)
        evaluator = ObjectiveEvaluator(
            prepared,
            objective,
            key_digits=work_dps + 8,
        )

        if previous_root_text is None:
            center = mp.mpf(center_text)
            half_width = mp.mpf(half_width_text)
        else:
            center = mp.mpf(previous_root_text)
            assert previous_dps is not None and previous_width_text is not None
            inherited_width = abs(mp.mpf(previous_width_text)) * 100
            precision_margin = mp.power(10, -max(6, previous_dps - 8))
            half_width = max(inherited_width, precision_margin)

        left, right, f_left, f_right, expansions = locate_minimum_bracket(
            evaluator,
            center,
            half_width,
            max_expansions=max_expansions,
        )
        x_tolerance = mp.power(10, -requested_dps)
        root_result = ridder_root(
            evaluator.stationary,
            left,
            right,
            f_left,
            f_right,
            x_tolerance=x_tolerance,
            max_steps=max_steps,
            evaluation_count=lambda: len(evaluator.cache),
        )

        # Always publish the canonical score, even when the faster raw-resultant
        # objective was selected for the location step.
        final_record = prepared.evaluate_score(root_result.root, second=True)

        stable_significant = None
        stable_decimals = None
        if previous_root_text is not None:
            stable_significant = significant_agreement(
                root_result.root, mp.mpf(previous_root_text), requested_dps
            )
            stable_decimals = decimal_place_agreement(
                root_result.root, mp.mpf(previous_root_text), requested_dps
            )

        payload = stage_payload(
            geometry=geometry,
            requested_dps=requested_dps,
            work_dps=work_dps,
            objective=objective,
            center=center,
            initial_half_width=half_width,
            root_result=root_result,
            bracket_expansions=expansions,
            final_record=final_record,
            elapsed_seconds=time.perf_counter() - started,
            stable_significant_digits=stable_significant,
            stable_decimal_places=stable_decimals,
        )
        root_text = mp_text(root_result.root, requested_dps + guard_digits - 2)
        width_text = mp_text(
            root_result.estimated_error, requested_dps + guard_digits - 2
        )
        return payload, root_text, width_text


def run_arithmetic_ladder(
    *,
    geometry: CameraGeometry,
    center_text: str,
    half_width_text: str,
    dps_values: Sequence[int],
    guard_digits: int,
    objective: str,
    max_steps: int,
    max_expansions: int,
) -> list[dict[str, Any]]:
    stages: list[dict[str, Any]] = []
    previous_root: str | None = None
    previous_width: str | None = None
    previous_dps: int | None = None
    for dps in dps_values:
        stage, previous_root, previous_width = run_precision_stage(
            geometry=geometry,
            center_text=center_text,
            half_width_text=half_width_text,
            requested_dps=dps,
            guard_digits=guard_digits,
            objective=objective,
            max_steps=max_steps,
            max_expansions=max_expansions,
            previous_root_text=previous_root,
            previous_dps=previous_dps,
            previous_width_text=previous_width,
        )
        stages.append(stage)
        previous_dps = dps
        print(
            f"  dps={dps:<4d} t={stage['root']}  "
            f"score={stage['score']}  evals={stage['function_evaluations']}  "
            f"elapsed={stage['elapsed_seconds']:.2f}s",
            flush=True,
        )
        agreement = stage["agreement_with_previous_stage"]
        if agreement["significant_digits"] is not None:
            print(
                "           arithmetic agreement: "
                f"{agreement['significant_digits']} significant digits, "
                f"{agreement['decimal_places']} decimal places",
                flush=True,
            )
    return stages


def run_cutoff_ladder(
    *,
    camera: int,
    cutoffs: Sequence[int],
    center_text: str,
    half_width_text: str,
    dps: int,
    guard_digits: int,
    objective: str,
    max_steps: int,
    max_expansions: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    previous_root: str | None = None
    for cutoff in cutoffs:
        geometry = build_geometry(camera, cutoff)
        # A new cutoff is a different finite function; retain the original
        # coarse half-width and let the bracketing routine expand if needed.
        stage, root_text, _ = run_precision_stage(
            geometry=geometry,
            center_text=previous_root or center_text,
            half_width_text=half_width_text,
            requested_dps=dps,
            guard_digits=guard_digits,
            objective=objective,
            max_steps=max_steps,
            max_expansions=max_expansions,
            previous_root_text=None,
            previous_dps=None,
            previous_width_text=None,
        )
        drift = None
        stable_significant = None
        stable_decimals = None
        if previous_root is not None:
            with mp.workdps(dps + guard_digits):
                drift_value = abs(mp.mpf(root_text) - mp.mpf(previous_root))
                drift = scientific_text(drift_value, 18)
                stable_significant = significant_agreement(
                    mp.mpf(root_text), mp.mpf(previous_root), dps
                )
                stable_decimals = decimal_place_agreement(
                    mp.mpf(root_text), mp.mpf(previous_root), dps
                )
        row = {
            "cutoff": cutoff,
            "root": stage["root"],
            "score": stage["score"],
            "resultant_norm": stage["resultant_norm"],
            "drift_from_previous_cutoff": drift,
            "agreement_with_previous_cutoff": {
                "significant_digits": stable_significant,
                "decimal_places": stable_decimals,
            },
            "elapsed_seconds": stage["elapsed_seconds"],
        }
        rows.append(row)
        previous_root = root_text
        print(
            f"  M={cutoff:<8d} t={row['root']}  score={row['score']}  "
            f"drift={drift or '-'}  elapsed={row['elapsed_seconds']:.2f}s",
            flush=True,
        )
    return rows


def self_test() -> int:
    print("SELF-TEST NATIVE-CARRY PRECISION LADDER")
    with mp.workdps(70):
        for camera in (2, 3, 4, 5):
            geometry = build_geometry(camera, 9)
            prepared = PreparedCamera(geometry)
            packed = prepared.evaluate_resultant(mp.mpf("14.125"))
            literal_x, literal_y = direct_resultant(prepared, mp.mpf("14.125"))
            tolerance = mp.mpf("1e-60")
            if abs(packed.resultant_x - literal_x) > tolerance:
                raise AssertionError(f"sparse/direct x mismatch for camera {camera}")
            if abs(packed.resultant_y - literal_y) > tolerance:
                raise AssertionError(f"sparse/direct y mismatch for camera {camera}")
            scored = prepared.evaluate_score(mp.mpf("14.125"))
            if scored.energy is None or scored.energy <= 0:
                raise AssertionError(f"nonpositive energy for camera {camera}")
            if scored.score is None or scored.score < 0:
                raise AssertionError(f"invalid score for camera {camera}")
        print("  PASS exact coefficient regrouping equals literal seed+bracket sum")
        print("  PASS C2 and saturated odd/even camera geometries")
        print("  PASS positive finite energy and nonnegative score")

        toy_calls = 0

        def toy(x: Any) -> Any:
            nonlocal toy_calls
            toy_calls += 1
            return x * x - 2

        left = mp.mpf("1")
        right = mp.mpf("2")
        result = ridder_root(
            toy,
            left,
            right,
            toy(left),
            toy(right),
            x_tolerance=mp.mpf("1e-55"),
            max_steps=100,
            evaluation_count=lambda: toy_calls,
        )
        if abs(result.root - mp.sqrt(2)) > mp.mpf("1e-50"):
            raise AssertionError("Ridder root solver failed")
        print("  PASS bracketed Ridder solver (no Newton)")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "High-precision finite native-carry R^2 minimum laboratory: "
            "float64/CUDA discovery plus mpmath Ridder refinement, without Newton"
        ),
    )
    parser.add_argument("--camera", type=int, default=3)
    parser.add_argument("--cutoff", type=int, default=16364)
    parser.add_argument(
        "--t",
        help="coarse candidate as an exact decimal string; omit only with --coarse",
    )
    parser.add_argument(
        "--half-width",
        default="0.0001",
        help="initial half-width around the coarse candidate",
    )
    parser.add_argument(
        "--dps",
        default=",".join(str(value) for value in DEFAULT_DPS),
        help="arithmetic precision ladder in decimal digits",
    )
    parser.add_argument("--guard-digits", type=int, default=20)
    parser.add_argument(
        "--objective",
        choices=("score", "resultant"),
        default="resultant",
        help=(
            "score matches the scanner exactly; resultant is faster and has the same exact zeros"
        ),
    )
    parser.add_argument("--max-steps", type=int, default=100)
    parser.add_argument("--max-expansions", type=int, default=12)
    parser.add_argument(
        "--cutoff-ladder",
        default="",
        help="optional finite-cutoff convergence ladder, e.g. 2048,4096,8192,16364",
    )
    parser.add_argument(
        "--cutoff-dps",
        type=int,
        default=60,
        help="precision used by the optional cutoff ladder",
    )
    parser.add_argument("--json-out", type=Path)
    parser.add_argument("--self-test", action="store_true")

    coarse = parser.add_argument_group("optional float64/CUDA coarse discovery")
    coarse.add_argument("--coarse", action="store_true")
    coarse.add_argument(
        "--operator",
        default=str(Path(__file__).with_name(DEFAULT_OPERATOR_NAME)),
    )
    coarse.add_argument("--coarse-tmin", default="92.48")
    coarse.add_argument("--coarse-tmax", default="92.50")
    coarse.add_argument("--coarse-grid", default="0.00001")
    coarse.add_argument(
        "--coarse-backend", choices=("auto", "cpu", "cuda"), default="auto"
    )
    coarse.add_argument("--workers", type=int, default=0)
    coarse.add_argument("--cpu-batch", type=int, default=16)
    coarse.add_argument("--state-block", type=int, default=32768)
    coarse.add_argument("--gpu-batch", type=int, default=8192)
    coarse.add_argument("--gpu-threads", type=int, default=256)
    coarse.add_argument("--prominence", type=float, default=0.0)
    coarse.add_argument("--prominence-window", type=int, default=16)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.self_test:
        return self_test()

    try:
        dps_values = parse_int_list(args.dps, name="dps", minimum=20)
        if tuple(sorted(set(dps_values))) != dps_values:
            raise ValueError("dps must be strictly increasing and distinct")
        cutoff_values = (
            parse_int_list(args.cutoff_ladder, name="cutoff-ladder", minimum=1)
            if args.cutoff_ladder.strip()
            else ()
        )
        if cutoff_values and tuple(sorted(set(cutoff_values))) != cutoff_values:
            raise ValueError("cutoff-ladder must be strictly increasing and distinct")
        if args.camera < 2:
            raise ValueError("camera must be >= 2")
        if args.cutoff < 1:
            raise ValueError("cutoff must be >= 1")
        if args.guard_digits < 8:
            raise ValueError("guard-digits must be >= 8")
        if args.max_steps < 2:
            raise ValueError("max-steps must be >= 2")
        if args.max_expansions < 0:
            raise ValueError("max-expansions must be >= 0")
        if args.cutoff_dps < 20:
            raise ValueError("cutoff-dps must be >= 20")
        if mp.mpf(args.half_width) <= 0:
            raise ValueError("half-width must be positive")
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    coarse_payload = None
    if args.coarse:
        candidate = discover_coarse_candidate(args)
        center_text = candidate.t_decimal
        half_width_text = candidate.grid_step
        coarse_payload = asdict(candidate)
        print("=" * 112)
        print(" FLOAT64/CUDA COARSE DISCOVERY")
        print("=" * 112)
        print(
            f"  camera={args.camera} cutoff={args.cutoff} backend={candidate.backend} "
            f"t={candidate.t_decimal} score={candidate.score:.17e} "
            f"elapsed={candidate.elapsed_seconds:.3f}s"
        )
    else:
        if args.t is None:
            raise SystemExit("provide --t or enable --coarse")
        center_text = args.t
        half_width_text = args.half_width

    # Validate exact decimal input before the expensive preparation.
    try:
        with mp.workdps(max(dps_values[0], 30)):
            candidate_value = mp.mpf(center_text)
            width_value = mp.mpf(half_width_text)
            if not mp.isfinite(candidate_value) or not mp.isfinite(width_value):
                raise ValueError("candidate and half-width must be finite")
            if width_value <= 0:
                raise ValueError("half-width must be positive")
    except Exception as exc:
        raise SystemExit(f"invalid decimal candidate/window: {exc}") from exc

    geometry = build_geometry(args.camera, args.cutoff)
    print("=" * 112)
    print(" NATIVE-CARRY PRIMITIVE R^2 PRECISION LADDER — NO NEWTON, NO INTERVAL ARITHMETIC")
    print("=" * 112)
    print(f" camera              : {geometry.camera}")
    print(f" geometry            : {geometry.geometry}")
    print(f" cutoff M            : {geometry.cutoff:,} centers")
    print(f" brackets            : {len(geometry.brackets):,}")
    print(f" coordinates N       : {geometry.coordinate_count:,}")
    print(f" unique states       : {len(geometry.unique_n):,}")
    print(f" largest center      : {geometry.largest_center:,}")
    print(f" coarse candidate    : {center_text}")
    print(f" initial half-width  : {half_width_text}")
    print(f" objective           : {args.objective}")
    print(f" precision ladder    : {dps_values}")
    print(f" mpmath              : {mp.__version__}")
    print("-" * 112)
    print(" Arithmetic ladder:", flush=True)

    all_started = time.perf_counter()
    arithmetic = run_arithmetic_ladder(
        geometry=geometry,
        center_text=center_text,
        half_width_text=half_width_text,
        dps_values=dps_values,
        guard_digits=args.guard_digits,
        objective=args.objective,
        max_steps=args.max_steps,
        max_expansions=args.max_expansions,
    )

    cutoff_rows: list[dict[str, Any]] = []
    if cutoff_values:
        print("-" * 112)
        print(f" Cutoff ladder at {args.cutoff_dps} dps:", flush=True)
        cutoff_rows = run_cutoff_ladder(
            camera=args.camera,
            cutoffs=cutoff_values,
            center_text=center_text,
            half_width_text=half_width_text,
            dps=args.cutoff_dps,
            guard_digits=args.guard_digits,
            objective=args.objective,
            max_steps=args.max_steps,
            max_expansions=args.max_expansions,
        )

    total_elapsed = time.perf_counter() - all_started
    final_stage = arithmetic[-1]
    arithmetic_agreement = final_stage["agreement_with_previous_stage"]
    cutoff_agreement = (
        cutoff_rows[-1]["agreement_with_previous_cutoff"] if cutoff_rows else None
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "status": "FINITE_NATIVE_CARRY_R2_PRECISION_AUDIT",
        "interpretation": {
            "newton": False,
            "interval_arithmetic": False,
            "root_solver": "bracketed Ridder on the stationary numerator",
            "operator_changed": False,
            "finite_cutoff": True,
            "tiny_score_is_not_a_digit_count": True,
        },
        "operator": {
            "state": "n^(-1/2) * (cos(-t log n), sin(-t log n))",
            "camera": geometry.camera,
            "geometry": geometry.geometry,
            "cutoff": geometry.cutoff,
            "half_range": geometry.half_range,
            "seed_count": len(geometry.seeds),
            "bracket_count": len(geometry.brackets),
            "coordinate_count": geometry.coordinate_count,
            "unique_state_count": len(geometry.unique_n),
            "largest_center": geometry.largest_center,
            "score": "norm(R)^2 / (N * sum_e norm(z_e)^2)",
            "post_bracket_map": "none",
            "coordinate_field": "R^2",
        },
        "coarse_candidate": coarse_payload
        or {
            "t_decimal": center_text,
            "half_width": half_width_text,
            "source": "explicit --t",
        },
        "arithmetic_ladder": arithmetic,
        "cutoff_ladder": cutoff_rows,
        "verdict": {
            "finite_cutoff_minimum": final_stage["root"],
            "finite_cutoff_score": final_stage["score"],
            "arithmetic_agreement_significant_digits": arithmetic_agreement[
                "significant_digits"
            ],
            "arithmetic_agreement_decimal_places": arithmetic_agreement[
                "decimal_places"
            ],
            "cutoff_agreement_significant_digits": (
                cutoff_agreement["significant_digits"]
                if cutoff_agreement is not None
                else None
            ),
            "cutoff_agreement_decimal_places": (
                cutoff_agreement["decimal_places"]
                if cutoff_agreement is not None
                else None
            ),
            "float64_ulp_at_root": final_stage["float64_projection"]["ulp_at_root"],
            "limiting_diagnostic": (
                "compare cutoff agreement against arithmetic agreement"
                if cutoff_rows
                else "run --cutoff-ladder to separate arithmetic digits from cutoff digits"
            ),
        },
        "runtime": {
            "elapsed_seconds": total_elapsed,
            "python": platform.python_version(),
            "platform": platform.platform(),
            "mpmath": mp.__version__,
            "script_sha256": file_sha256(Path(__file__).resolve()),
        },
    }

    print("-" * 112)
    print(" VERDICT")
    print(f" finite-cutoff minimum : {final_stage['root']}")
    print(f" score                  : {final_stage['score']}")
    print(
        " float64 ULP at t       : "
        f"{final_stage['float64_projection']['ulp_at_root']}"
    )
    if arithmetic_agreement["significant_digits"] is not None:
        print(
            " arithmetic agreement  : "
            f"{arithmetic_agreement['significant_digits']} significant digits "
            f"({arithmetic_agreement['decimal_places']} decimal places)"
        )
    if cutoff_agreement is not None:
        print(
            " cutoff agreement      : "
            f"{cutoff_agreement['significant_digits']} significant digits "
            f"({cutoff_agreement['decimal_places']} decimal places)"
        )
        print(" practical precision    : the smaller of arithmetic and cutoff agreement")
    else:
        print(" cutoff agreement      : not measured; add --cutoff-ladder")
    print(f" total elapsed          : {total_elapsed:.2f}s")

    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.json_out.with_suffix(args.json_out.suffix + ".tmp")
        temporary.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(args.json_out)
        print(f" JSON written           : {args.json_out}")
    print("=" * 112)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
