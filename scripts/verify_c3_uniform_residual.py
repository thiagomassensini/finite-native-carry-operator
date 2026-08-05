#!/usr/bin/env python3
"""Verify the decomposed cutoff-uniform route for the C3 residual."""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any, Sequence


SCHEMA = "org.native-carry.c3-uniform-residual-decomposition/v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _integer_bounds(enclosure: dict[str, Any]) -> tuple[int, int, int]:
    midpoint = int(enclosure["midpoint_integer"])
    radius = int(enclosure["radius_integer"])
    exponent = int(enclosure["exponent10"])
    if radius < 0:
        raise AssertionError("an enclosure radius cannot be negative")
    return midpoint - radius, midpoint + radius, exponent


def _verify_upper_record(record: dict[str, Any]) -> None:
    _, enclosure_upper, exponent = _integer_bounds(record["formula_enclosure"])
    upper = record["certified_upper"]
    if int(upper["significand_integer"]) != enclosure_upper:
        raise AssertionError("certified upper is not its enclosure endpoint")
    if int(upper["exponent10"]) != exponent:
        raise AssertionError("certified upper exponent mismatch")
    if enclosure_upper <= 0:
        raise AssertionError("certified upper must be positive")


def _decimal_upper(record: dict[str, Any]) -> Decimal:
    upper = record["certified_upper"]
    return Decimal(upper["significand_integer"]) * (
        Decimal(10) ** int(upper["exponent10"])
    )


def verify_uniform_residual_decomposition(
    path: Path,
    source_root: Path | None = None,
    *,
    recompute: bool = False,
) -> dict[str, Any]:
    getcontext().prec = 1000
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload["schema"] != SCHEMA:
        raise AssertionError("unsupported uniform-residual schema")
    scope = payload["scope"]
    if scope["coordinate_field"] != "R^2":
        raise AssertionError("the residual decomposition is not real-plane native")
    if scope["theory_scope"] != "native_real_operator_only":
        raise AssertionError("the decomposition left the native real theory")
    if int(scope["camera"]) != 3 or scope["shared_stationary_point"] is not True:
        raise AssertionError("the decomposition does not fix the C3 stationary point")

    anchor = payload["anchor"]
    _verify_upper_record(anchor["common_limiting_velocity_norm_upper_V"])
    _verify_upper_record(anchor["space_factor_P6_T_upper"])
    if Decimal(anchor["domain_lower"]) >= Decimal(anchor["domain_upper"]):
        raise AssertionError("the anchor domain is empty")

    entries = payload["entries"]
    if len(entries) < 2:
        raise AssertionError("the finite decomposition needs at least two entries")
    cutoffs = [int(entry["cutoff"]) for entry in entries]
    if cutoffs != sorted(set(cutoffs)) or cutoffs[0] < 2:
        raise AssertionError("the decomposition cutoffs are invalid")

    record_keys = (
        "core_residual_Q_M_upper",
        "sharp_oriented_tail_eta_M_upper",
        "polynomial_tail_witness_upper",
        "common_velocity_times_radius_upper",
        "decomposed_stationary_norm_upper",
        "decomposed_stationary_energy_upper",
        "original_stationary_norm_upper",
    )
    root = source_root or path.resolve().parents[2]
    core_values: list[Decimal] = []
    radii: list[Decimal] = []
    for entry in entries:
        for key in record_keys:
            _verify_upper_record(entry[key])
        certificate_path = root / entry["oriented_certificate_path"]
        if _sha256(certificate_path) != entry["oriented_certificate_sha256"]:
            raise AssertionError(f"oriented certificate digest mismatch: {certificate_path}")

        core = _decimal_upper(entry["core_residual_Q_M_upper"])
        tail = _decimal_upper(entry["sharp_oriented_tail_eta_M_upper"])
        polynomial = _decimal_upper(entry["polynomial_tail_witness_upper"])
        localization = _decimal_upper(
            entry["common_velocity_times_radius_upper"]
        )
        norm = _decimal_upper(entry["decomposed_stationary_norm_upper"])
        energy = _decimal_upper(entry["decomposed_stationary_energy_upper"])
        original = _decimal_upper(entry["original_stationary_norm_upper"])
        if tail > polynomial:
            raise AssertionError("the polynomial witness does not bound the tail")
        if norm < core + tail + localization:
            raise AssertionError("the decomposed norm does not cover its components")
        if energy < norm * norm:
            raise AssertionError("the energy envelope does not cover the squared norm")
        if entry["decomposition_improves_original_norm_upper"] != (
            norm < original
        ):
            raise AssertionError("the original-bound comparison is inconsistent")
        core_values.append(core)
        radii.append(Decimal(entry["radius"]))

    core_contracts = all(
        current < previous
        for previous, current in zip(core_values, core_values[1:])
    )
    radii_contract = all(
        current < previous for previous, current in zip(radii, radii[1:])
    )
    expected_status = {
        "sharp_tail_remainder_bounded_by_polynomial_witness": True,
        "tail_component_tends_to_zero_for_fixed_T": True,
        "localization_component_tends_to_zero_if_radii_tend_to_zero": True,
        "core_residual_Q_M_strictly_contracts_on_finite_entries": core_contracts,
        "radii_strictly_contract_on_finite_entries": radii_contract,
        "core_residual_Q_M_infinite_vanishing_family_certified": False,
        "radius_infinite_vanishing_family_certified": False,
        "limiting_vector_zero_certified": False,
    }
    if payload["component_status"] != expected_status:
        raise AssertionError("the component claims exceed the finite ledger")
    if not (core_contracts and radii_contract):
        raise AssertionError("the finite core or radius sequence does not contract")

    decomposition = payload["decomposition"]
    if decomposition["norm_bound"] != (
        "||R_infinity(t_*)|| <= Q_M + eta_M + V*r_M"
    ):
        raise AssertionError("unexpected residual decomposition")
    if decomposition["polynomial_witness"] != (
        "eta_M <= P6(T)*(16/1485)/(M+1)^5 for every integer M >= 2"
    ):
        raise AssertionError("unexpected polynomial tail witness")
    lean_bridge = payload["lean_bridge"]
    if lean_bridge["contract"] != (
        "DecomposedVanishingLimitResidualCertificate.witness_zero"
    ):
        raise AssertionError("unexpected decomposed Lean bridge")
    if lean_bridge["tail_limit_theorem"] != (
        "tendsto_zero_of_le_polynomialTailEnvelope"
    ):
        raise AssertionError("unexpected polynomial-tail Lean theorem")

    for source in payload["source_hashes"].values():
        source_path = root / source["path"]
        if _sha256(source_path) != source["sha256"]:
            raise AssertionError(f"source digest mismatch: {source_path}")

    if recompute:
        import sys

        root_text = str(root)
        if root_text not in sys.path:
            sys.path.insert(0, root_text)
        from certification.c3_uniform_residual import (
            build_uniform_residual_decomposition,
        )

        rebuilt = build_uniform_residual_decomposition(
            [root / entry["oriented_certificate_path"] for entry in entries],
            source_root=root,
        )
        if rebuilt != payload:
            differing = sorted(
                key for key in rebuilt if rebuilt.get(key) != payload.get(key)
            )
            raise AssertionError(
                f"recomputed uniform-residual mismatch: {differing}"
            )

    return payload


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args(argv)
    payload = verify_uniform_residual_decomposition(
        args.ledger, recompute=args.recompute
    )
    print(f"verified={args.ledger}")
    print(f"recomputed={str(args.recompute).lower()}")
    for key, value in payload["component_status"].items():
        print(f"{key}={str(value).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
