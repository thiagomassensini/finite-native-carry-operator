from __future__ import annotations

from decimal import Decimal, getcontext
from pathlib import Path

from scripts.verify_c3_stationary_localization import (
    verify_stationary_localization,
)


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "certificates/c3/c3_stationary_localization.json"


def _upper(record: dict[str, object]) -> Decimal:
    certified = record["certified_upper"]
    assert isinstance(certified, dict)
    return Decimal(str(certified["significand_integer"])) * (
        Decimal(10) ** int(str(certified["exponent10"]))
    )


def test_stationary_localization_recomputes_exactly() -> None:
    payload = verify_stationary_localization(LEDGER, ROOT, recompute=True)
    assert [entry["cutoff"] for entry in payload["entries"]] == [
        16364,
        65536,
        131072,
    ]
    claims = payload["claims"]
    assert claims["finite_stationary_localization_radii_certified"]
    assert claims["ideal_corrected_root_radius_witness_tends_to_zero"]
    assert not claims[
        "corrected_stationary_root_family_constructed_for_all_cutoffs"
    ]
    assert not claims["limiting_vector_zero_certified"]


def test_derived_radii_fit_inside_every_oriented_interval() -> None:
    getcontext().prec = 1000
    payload = verify_stationary_localization(LEDGER, ROOT)
    for entry in payload["entries"]:
        derived = _upper(entry["derived_localization_radius_upper"])
        certified = Decimal(entry["certified_radius"])
        assert derived <= certified
        assert entry["derived_radius_within_certified_radius"]
        assert _upper(
            entry["resultant_remainder_polynomial_witness_upper"]
        ) >= _upper(entry["sharp_resultant_remainder_eta0_upper"])
        assert _upper(
            entry["velocity_remainder_polynomial_witness_upper"]
        ) >= _upper(entry["sharp_velocity_remainder_eta1_upper"])
