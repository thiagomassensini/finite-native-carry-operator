from __future__ import annotations

from pathlib import Path

from scripts.verify_c3_contraction_ladder import verify_contraction_ladder


ROOT = Path(__file__).resolve().parents[1]
LADDER = ROOT / "certificates/c3/c3_oriented_contraction_ladder.json"


def test_contraction_ladder_recomputes_exactly() -> None:
    payload = verify_contraction_ladder(LADDER, ROOT, recompute=True)
    assert [entry["cutoff"] for entry in payload["entries"]] == [
        16364,
        65536,
        131072,
    ]
    assert payload["claims"][
        "all_entries_bound_the_same_limiting_stationary_point"
    ]
    assert payload["claims"][
        "resultant_energy_upper_bounds_strictly_contract_across_entries"
    ]
    assert not payload["claims"]["infinite_vanishing_bound_family_certified"]
    assert not payload["claims"]["limiting_vector_zero_certified"]


def test_energy_uppers_are_exact_squares_of_norm_uppers() -> None:
    payload = verify_contraction_ladder(LADDER, ROOT)
    for entry in payload["entries"]:
        norm = entry["resultant_norm_upper_at_shared_stationary_point"]
        energy = entry["resultant_energy_upper_at_shared_stationary_point"]
        assert int(energy["significand_integer"]) == int(
            norm["significand_integer"]
        ) ** 2
        assert int(energy["exponent10"]) == 2 * int(norm["exponent10"])
