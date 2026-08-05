#!/usr/bin/env python3
"""Decompose C3 stationary residual bounds into three uniform obligations.

For a certified limiting stationary point t_* in a nested interval
[c_M-r_M,c_M+r_M], the oriented tail theorem gives

    ||R_inf(t_*)|| <= Q_M + eta_M + V * r_M.

Here Q_M is the norm of the finite sum plus the fifth-order boundary jet at
c_M, eta_M is the oriented sixth-derivative remainder, and V is one velocity
cap certified on the largest interval.  The resultant remainder has the
closed form

    eta_M(T) = P_6(T) * (
        7/660 * (3(M+1))^(-11/2)
        + 1/5940 * (3M-1)^(-11/2)
    ),

and, for M >= 2, is bounded by

    P_6(T) * 16/1485 / (M+1)^5.

Thus the tail term has a cutoff-uniform proof of convergence to zero.  The
ledger deliberately leaves convergence of Q_M as a separate obligation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal
from pathlib import Path
from typing import Any, Sequence

import flint
from flint import arb, ctx

from certification.c3_oriented_tail import (
    _space_factor,
    evaluate_oriented_c3_limit,
    oriented_c3_tail_enclosure,
)
from certification.c3_tail import _plane_norm_upper, _positive_upper_record, _upper_ball
from certification.real_interval import (
    build_sparse_geometry,
    evaluate_real_operator,
    prepare_terms,
)
from scripts.verify_c3_oriented_tail_certificate import verify_oriented_certificate


SCHEMA = "org.native-carry.c3-uniform-residual-decomposition/v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _decimal_upper(record: dict[str, Any]) -> Decimal:
    upper = record["certified_upper"]
    return Decimal(upper["significand_integer"]) * (
        Decimal(10) ** int(upper["exponent10"])
    )


def build_uniform_residual_decomposition(
    certificate_paths: Sequence[Path],
    *,
    source_root: Path | None = None,
) -> dict[str, Any]:
    if len(certificate_paths) < 2:
        raise ValueError("the decomposition needs at least two certificates")
    root = source_root or Path(__file__).resolve().parents[1]

    loaded: list[tuple[Path, dict[str, Any]]] = []
    for path in certificate_paths:
        resolved = path if path.is_absolute() else root / path
        certificate = verify_oriented_certificate(resolved, root, recompute=False)
        loaded.append((resolved, certificate))
    loaded.sort(key=lambda item: int(item[1]["operator"]["finite_cutoff"]))

    limit_objects = {certificate["scope"]["limit_object"] for _, certificate in loaded}
    if len(limit_objects) != 1:
        raise RuntimeError("the certificates do not enclose the same C3 limit")

    anchor_path, anchor = loaded[0]
    anchor_lower = Decimal(anchor["domain"]["exact_decimal_lower"])
    anchor_upper = Decimal(anchor["domain"]["exact_decimal_upper"])
    previous_lower = anchor_lower
    previous_upper = anchor_upper
    for _, certificate in loaded[1:]:
        lower = Decimal(certificate["domain"]["exact_decimal_lower"])
        upper = Decimal(certificate["domain"]["exact_decimal_upper"])
        if not (previous_lower <= lower <= upper <= previous_upper):
            raise RuntimeError("the oriented intervals are not nested")
        previous_lower, previous_upper = lower, upper

    precision = max(
        int(certificate["arithmetic"]["decimal_digits"])
        for _, certificate in loaded
    )
    ctx.dps = precision
    digits = precision + 12
    time_cap_text = anchor["domain"]["time_abs_upper"]
    time_cap = arb(time_cap_text)

    anchor_cutoff = int(anchor["operator"]["finite_cutoff"])
    anchor_geometry = build_sparse_geometry(3, anchor_cutoff)
    anchor_prepared = prepare_terms(anchor_geometry)
    anchor_domain = arb(
        arb(anchor["domain"]["requested_center"]),
        arb(anchor["domain"]["requested_radius"]),
    )
    anchor_evaluation = evaluate_oriented_c3_limit(
        time=anchor_domain,
        prepared=anchor_prepared,
        cutoff=anchor_cutoff,
        time_abs_upper=time_cap,
    )
    velocity_cap_record = _plane_norm_upper(
        anchor_evaluation.derivative_x,
        anchor_evaluation.derivative_y,
        digits,
    )
    velocity_cap = _upper_ball(velocity_cap_record)
    space_factor_six_record = _positive_upper_record(
        _space_factor(6, time_cap), digits
    )
    space_factor_six = _upper_ball(space_factor_six_record)

    entries: list[dict[str, Any]] = []
    for path, certificate in loaded:
        cutoff = int(certificate["operator"]["finite_cutoff"])
        center_text = certificate["domain"]["requested_center"]
        radius_text = certificate["domain"]["requested_radius"]
        radius = arb(radius_text)
        geometry = build_sparse_geometry(3, cutoff)
        prepared = prepare_terms(geometry)

        finite_center = evaluate_real_operator(
            arb(center_text), prepared, second=False
        )
        tail_center = oriented_c3_tail_enclosure(
            cutoff=cutoff,
            time=arb(center_text),
            time_abs_upper=time_cap,
            time_order=0,
        )
        corrected_x = finite_center.resultant_x + tail_center.approximation_x
        corrected_y = finite_center.resultant_y + tail_center.approximation_y
        core_record = _plane_norm_upper(corrected_x, corrected_y, digits)
        sharp_tail_record = _positive_upper_record(
            tail_center.remainder_norm, digits
        )
        polynomial_tail_record = _positive_upper_record(
            space_factor_six
            * arb(16)
            / 1485
            / (arb(cutoff + 1) ** 5),
            digits,
        )
        if _decimal_upper(sharp_tail_record) > _decimal_upper(
            polynomial_tail_record
        ):
            raise RuntimeError("the polynomial tail witness did not dominate")

        localization_record = _positive_upper_record(
            velocity_cap * radius, digits
        )
        norm_envelope_record = _positive_upper_record(
            _upper_ball(core_record)
            + _upper_ball(sharp_tail_record)
            + _upper_ball(localization_record),
            digits,
        )
        energy_envelope_record = _positive_upper_record(
            _upper_ball(norm_envelope_record) ** 2,
            digits,
        )
        original_record = certificate["stationary_point_bounds"][
            "resultant_norm_upper"
        ]

        entries.append(
            {
                "cutoff": cutoff,
                "oriented_certificate_path": str(path.relative_to(root)),
                "oriented_certificate_sha256": _sha256(path),
                "center": center_text,
                "radius": radius_text,
                "corrected_center_resultant_x_display": corrected_x.str(
                    digits, radius=True, more=True
                ),
                "corrected_center_resultant_y_display": corrected_y.str(
                    digits, radius=True, more=True
                ),
                "core_residual_Q_M_upper": core_record,
                "sharp_oriented_tail_eta_M_upper": sharp_tail_record,
                "polynomial_tail_witness_upper": polynomial_tail_record,
                "common_velocity_times_radius_upper": localization_record,
                "decomposed_stationary_norm_upper": norm_envelope_record,
                "decomposed_stationary_energy_upper": energy_envelope_record,
                "original_stationary_norm_upper": original_record,
                "decomposition_improves_original_norm_upper": (
                    _decimal_upper(norm_envelope_record)
                    < _decimal_upper(original_record)
                ),
            }
        )

    core_values = [
        _decimal_upper(entry["core_residual_Q_M_upper"]) for entry in entries
    ]
    radii = [Decimal(entry["radius"]) for entry in entries]
    core_strictly_contracts = all(
        current < previous
        for previous, current in zip(core_values, core_values[1:])
    )
    radii_strictly_contract = all(
        current < previous for previous, current in zip(radii, radii[1:])
    )
    if not (core_strictly_contracts and radii_strictly_contract):
        raise RuntimeError("the finite residual components do not contract")

    source_paths = {
        "decomposition_builder": Path(__file__).resolve(),
        "oriented_tail_certifier": root / "certification/c3_oriented_tail.py",
        "finite_interval_evaluator": root / "certification/real_interval.py",
        "lean_contract": root
        / "FiniteNativeCarryOperator/Certification/Contract.lean",
    }
    return {
        "schema": SCHEMA,
        "status": "CERTIFIED_FINITE_C3_UNIFORM_RESIDUAL_DECOMPOSITION",
        "scope": {
            "coordinate_field": "R^2",
            "theory_scope": "native_real_operator_only",
            "camera": 3,
            "limit_object": next(iter(limit_objects)),
            "shared_stationary_point": True,
        },
        "anchor": {
            "oriented_certificate_path": str(anchor_path.relative_to(root)),
            "domain_lower": anchor["domain"]["exact_decimal_lower"],
            "domain_upper": anchor["domain"]["exact_decimal_upper"],
            "time_abs_upper_T": time_cap_text,
            "common_limiting_velocity_norm_upper_V": velocity_cap_record,
            "space_factor_P6_T_upper": space_factor_six_record,
        },
        "decomposition": {
            "norm_bound": "||R_infinity(t_*)|| <= Q_M + eta_M + V*r_M",
            "Q_M": "norm of the finite C3 sum plus the fifth-order oriented boundary jet at the interval center",
            "eta_M": "sixth-space-derivative oriented-tail remainder at the common time cap T",
            "V_r_M": "one common limiting-velocity cap on the anchor interval times the entry radius",
            "sharp_eta_formula": "P6(T)*(7/660*(3*(M+1))^(-11/2) + 1/5940*(3*M-1)^(-11/2))",
            "polynomial_witness": "eta_M <= P6(T)*(16/1485)/(M+1)^5 for every integer M >= 2",
        },
        "entries": entries,
        "component_status": {
            "sharp_tail_remainder_bounded_by_polynomial_witness": True,
            "tail_component_tends_to_zero_for_fixed_T": True,
            "localization_component_tends_to_zero_if_radii_tend_to_zero": True,
            "core_residual_Q_M_strictly_contracts_on_finite_entries": (
                core_strictly_contracts
            ),
            "radii_strictly_contract_on_finite_entries": radii_strictly_contract,
            "core_residual_Q_M_infinite_vanishing_family_certified": False,
            "radius_infinite_vanishing_family_certified": False,
            "limiting_vector_zero_certified": False,
        },
        "remaining_obligation": {
            "primary": "construct a cutoff-uniform family of corrected centers and prove Q_M tends to zero",
            "secondary": "prove the certified nested radii form a cutoff-uniform sequence tending to zero",
            "completed": "the analytic tail component eta_M has an explicit polynomial vanishing witness",
        },
        "lean_bridge": {
            "contract": "DecomposedVanishingLimitResidualCertificate.witness_zero",
            "tail_limit_theorem": "tendsto_zero_of_le_polynomialTailEnvelope",
            "required_component_limits": [
                "Q_M tends to zero",
                "eta_M tends to zero",
                "V*r_M tends to zero",
            ],
        },
        "arithmetic": {
            "backend": "python-flint Arb real balls",
            "python_flint_version": flint.__version__,
            "flint_version": flint.__FLINT_VERSION__,
            "decimal_digits": precision,
            "binary_precision_bits": ctx.prec,
            "directed_rounding": True,
        },
        "source_hashes": {
            name: {
                "path": str(path.relative_to(root)),
                "sha256": _sha256(path),
            }
            for name, path in source_paths.items()
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the decomposed C3 stationary-residual ledger"
    )
    parser.add_argument("certificates", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    payload = build_uniform_residual_decomposition(args.certificates)
    root = Path(__file__).resolve().parents[1]
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"uniform residual decomposition written to {output}")
    for key, value in payload["component_status"].items():
        print(f"{key}={str(value).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
