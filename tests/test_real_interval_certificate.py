from __future__ import annotations

import json
from pathlib import Path

from certification.real_interval import build_sparse_geometry
from laboratory.native_carry_precision_ladder import build_geometry
from scripts.verify_interval_certificate import verify_certificate


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates/c3/c3_m16364_raw_minimum.json"


def test_sparse_geometry_matches_precision_laboratory() -> None:
    for camera in range(2, 10):
        for cutoff in (1, 17):
            certified = build_sparse_geometry(camera, cutoff)
            laboratory = build_geometry(camera, cutoff)
            assert certified.terms == tuple(
                zip(laboratory.unique_n, laboratory.coefficients, strict=True)
            )
            assert certified.coordinate_count == laboratory.coordinate_count
            assert certified.largest_center == laboratory.largest_center


def test_c2_is_the_radius_one_sector_of_c4() -> None:
    cutoff = 31
    c2 = build_sparse_geometry(2, cutoff)

    coefficients: dict[int, int] = {1: 1}
    for m in range(1, cutoff + 1):
        center = 4 * m
        coefficients[center - 1] = coefficients.get(center - 1, 0) + 1
        coefficients[center] = coefficients.get(center, 0) - 2
        coefficients[center + 1] = coefficients.get(center + 1, 0) + 1
    radius_one = tuple(sorted((n, c) for n, c in coefficients.items() if c))
    assert c2.terms == radius_one


def test_checked_certificate_is_self_consistent() -> None:
    certificate = verify_certificate(CERTIFICATE, ROOT, recompute=True)
    assert certificate["operator"]["camera"] == 3
    assert certificate["operator"]["cutoff"] == 16364
    assert certificate["claims"] == {
        "finite_vector_zero_absent_from_domain": True,
        "limiting_zero_certified": False,
        "strict_finite_minimum_in_domain": True,
        "unique_stationary_point_in_domain": True,
    }


def test_certificate_is_valid_json() -> None:
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    assert payload["status"] == "REAL_INTERVAL_CERTIFIED_FINITE_MINIMUM"
