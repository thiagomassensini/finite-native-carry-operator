from __future__ import annotations

from decimal import Decimal, getcontext
from pathlib import Path

from scripts.verify_c3_uniform_residual import (
    verify_uniform_residual_decomposition,
)


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "certificates/c3/c3_uniform_residual_decomposition.json"


def _upper(record: dict[str, object]) -> Decimal:
    certified = record["certified_upper"]
    assert isinstance(certified, dict)
    return Decimal(str(certified["significand_integer"])) * (
        Decimal(10) ** int(str(certified["exponent10"]))
    )


def test_uniform_residual_decomposition_recomputes_exactly() -> None:
    payload = verify_uniform_residual_decomposition(
        LEDGER, ROOT, recompute=True
    )
    assert [entry["cutoff"] for entry in payload["entries"]] == [
        16364,
        65536,
        131072,
    ]
    status = payload["component_status"]
    assert status["tail_component_tends_to_zero_for_fixed_T"]
    assert status["core_residual_Q_M_strictly_contracts_on_finite_entries"]
    assert not status[
        "core_residual_Q_M_infinite_vanishing_family_certified"
    ]
    assert not status["limiting_vector_zero_certified"]


def test_each_component_is_exposed_without_hiding_the_core_residual() -> None:
    getcontext().prec = 1000
    payload = verify_uniform_residual_decomposition(LEDGER, ROOT)
    for entry in payload["entries"]:
        core = _upper(entry["core_residual_Q_M_upper"])
        tail = _upper(entry["sharp_oriented_tail_eta_M_upper"])
        localization = _upper(entry["common_velocity_times_radius_upper"])
        norm = _upper(entry["decomposed_stationary_norm_upper"])
        assert core > 0
        assert norm >= core + tail + localization
        assert _upper(entry["polynomial_tail_witness_upper"]) >= tail
