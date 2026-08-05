from __future__ import annotations

from pathlib import Path

from flint import arb, ctx

from certification.c3_oriented_tail import (
    mixed_space_time_derivative,
    oriented_c3_tail_enclosure,
)
from certification.c3_tail import c3_tail_majorants
from certification.real_interval import (
    PreparedTerm,
    build_sparse_geometry,
    evaluate_real_operator,
)
from scripts.verify_c3_oriented_tail_certificate import verify_oriented_certificate


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates/c3/c3_m16364_oriented_limit_minimum.json"
REFINED_CERTIFICATE = ROOT / "certificates/c3/c3_m65536_oriented_limit_minimum.json"
THIRD_CERTIFICATE = ROOT / "certificates/c3/c3_m131072_oriented_limit_minimum.json"


def _partial_c3_tail(start_cutoff: int, end_cutoff: int) -> tuple[PreparedTerm, ...]:
    start = dict(build_sparse_geometry(3, start_cutoff).terms)
    end = dict(build_sparse_geometry(3, end_cutoff).terms)
    terms: list[PreparedTerm] = []
    for n in sorted(set(start) | set(end)):
        coefficient = end.get(n, 0) - start.get(n, 0)
        if coefficient:
            n_ball = arb(n)
            terms.append(
                PreparedTerm(
                    amplitude=arb(coefficient) / n_ball.sqrt(),
                    log_n=n_ball.log(),
                )
            )
    return tuple(terms)


def test_time_derivative_convention_matches_the_real_rotation() -> None:
    ctx.dps = 70
    x_value = arb(17)
    time = arb("3.25")
    logarithm = x_value.log()
    amplitude = 1 / x_value.sqrt()
    angle = -time * logarithm
    base = (amplitude * angle.cos(), amplitude * angle.sin())
    expected = (
        base,
        (amplitude * logarithm * angle.sin(), -amplitude * logarithm * angle.cos()),
        (-logarithm * logarithm * base[0], -logarithm * logarithm * base[1]),
    )
    for order in range(3):
        actual = mixed_space_time_derivative(
            x_value=x_value,
            time=time,
            space_order=0,
            time_order=order,
        )
        assert (actual[0] - expected[order][0]).contains(0)
        assert (actual[1] - expected[order][1]).contains(0)


def test_oriented_tail_agrees_with_an_independent_long_partial_sum() -> None:
    ctx.dps = 70
    cutoff = 20
    comparison_cutoff = 4000
    time = arb("2.75")
    partial = evaluate_real_operator(
        time,
        _partial_c3_tail(cutoff, comparison_cutoff),
        second=True,
    )
    assert partial.second_x is not None and partial.second_y is not None
    reference = (
        (partial.resultant_x, partial.resultant_y),
        (partial.derivative_x, partial.derivative_y),
        (partial.second_x, partial.second_y),
    )
    remaining = c3_tail_majorants(
        cutoff=comparison_cutoff, time_abs_upper=time
    )
    remaining_bounds = (
        remaining.resultant,
        remaining.first_time_derivative,
        remaining.second_time_derivative,
    )

    for order in range(3):
        oriented = oriented_c3_tail_enclosure(
            cutoff=cutoff,
            time=time,
            time_abs_upper=time,
            time_order=order,
        )
        independent_x = arb(reference[order][0], remaining_bounds[order])
        independent_y = arb(reference[order][1], remaining_bounds[order])
        assert (oriented.enclosure_x - independent_x).contains(0)
        assert (oriented.enclosure_y - independent_y).contains(0)


def test_oriented_certificate_recomputes_exactly() -> None:
    certificate = verify_oriented_certificate(CERTIFICATE, ROOT, recompute=True)
    assert certificate["operator"]["camera"] == 3
    assert certificate["domain"]["requested_radius"] == "5e-16"
    assert certificate["claims"]["unique_limiting_stationary_point_in_domain"]
    assert certificate["claims"]["stationary_resultant_has_certified_small_norm"]
    assert not certificate["claims"]["limiting_vector_zero_certified"]


def test_refined_oriented_certificate_recomputes_exactly() -> None:
    certificate = verify_oriented_certificate(
        REFINED_CERTIFICATE, ROOT, recompute=True
    )
    assert certificate["operator"]["finite_cutoff"] == 65536
    assert certificate["domain"]["requested_radius"] == "2e-20"
    assert certificate["claims"]["unique_limiting_stationary_point_in_domain"]
    assert certificate["claims"]["stationary_resultant_has_certified_small_norm"]
    assert not certificate["claims"]["limiting_vector_zero_certified"]


def test_third_oriented_certificate_recomputes_exactly() -> None:
    certificate = verify_oriented_certificate(
        THIRD_CERTIFICATE, ROOT, recompute=True
    )
    assert certificate["operator"]["finite_cutoff"] == 131072
    assert certificate["domain"]["requested_radius"] == "4e-22"
    assert certificate["claims"]["unique_limiting_stationary_point_in_domain"]
    assert certificate["claims"]["stationary_resultant_has_certified_small_norm"]
    assert not certificate["claims"]["limiting_vector_zero_certified"]
