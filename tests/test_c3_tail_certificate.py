from __future__ import annotations

from pathlib import Path

from flint import arb, ctx

from certification.c3_tail import c3_tail_majorants
from certification.real_interval import (
    PreparedTerm,
    build_sparse_geometry,
    evaluate_real_operator,
)
from scripts.verify_c3_tail_certificate import verify_c3_tail_certificate


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates/c3/c3_m16364_tail_limit_minimum.json"


def _partial_c3_tail(start_cutoff: int, end_cutoff: int) -> tuple[PreparedTerm, ...]:
    start = dict(build_sparse_geometry(3, start_cutoff).terms)
    end = dict(build_sparse_geometry(3, end_cutoff).terms)
    indices = sorted(set(start) | set(end))
    result: list[PreparedTerm] = []
    for n in indices:
        coefficient = end.get(n, 0) - start.get(n, 0)
        if coefficient:
            n_ball = arb(n)
            result.append(
                PreparedTerm(
                    amplitude=arb(coefficient) / n_ball.sqrt(),
                    log_n=n_ball.log(),
                )
            )
    return tuple(result)


def _point_norm(x_value: object, y_value: object) -> object:
    return (x_value * x_value + y_value * y_value).sqrt()


def test_explicit_tail_bounds_dominate_a_long_partial_tail() -> None:
    ctx.dps = 70
    cutoff = 20
    time = arb("92.4919")
    evaluation = evaluate_real_operator(
        time, _partial_c3_tail(cutoff, 4000), second=True
    )
    majorants = c3_tail_majorants(cutoff=cutoff, time_abs_upper=time)
    assert evaluation.second_x is not None and evaluation.second_y is not None

    assert _point_norm(evaluation.resultant_x, evaluation.resultant_y) < (
        majorants.resultant
    )
    assert _point_norm(evaluation.derivative_x, evaluation.derivative_y) < (
        majorants.first_time_derivative
    )
    assert _point_norm(evaluation.second_x, evaluation.second_y) < (
        majorants.second_time_derivative
    )


def test_tail_majorants_contract_when_the_cutoff_doubles() -> None:
    ctx.dps = 60
    first = c3_tail_majorants(cutoff=100, time_abs_upper=arb(93))
    second = c3_tail_majorants(cutoff=200, time_abs_upper=arb(93))
    assert second.resultant < first.resultant
    assert second.first_time_derivative < first.first_time_derivative
    assert second.second_time_derivative < first.second_time_derivative


def test_c3_tail_certificate_recomputes_exactly() -> None:
    certificate = verify_c3_tail_certificate(
        CERTIFICATE, ROOT, recompute=True
    )
    assert certificate["operator"]["camera"] == 3
    assert certificate["operator"]["finite_cutoff"] == 16364
    assert certificate["claims"]["unique_limiting_stationary_point_in_domain"]
    assert certificate["claims"]["strict_limiting_minimum_in_domain"]
    assert certificate["claims"]["limiting_first_derivative_nonzero_on_domain"]
    assert certificate["claims"][
        "vector_zero_reduced_to_determinant_at_stationary_point"
    ]
    assert not certificate["claims"]["limiting_vector_zero_certified"]
