#!/usr/bin/env python3
"""Industrial all-base carry Green/Bessel radial atlas (CPU sieve + CUDA).

This program separates three mathematical objects that were previously mixed:

1. Camera family
   The primary atlas is every positional base b >= 2.  Odd, even, and prime
   subatlases are accumulated in parallel for comparison.  The prime subatlas
   reproduces the index family used by the existing Lean Green/Bessel theorem;
   the all-base atlas is an explicit numerical/generalized finite ledger.

2. Radial Green readout

       a_b(delta) = b^(-1/2) * (b^delta - b^(-delta))

   and therefore

       readout_energy = sum_b a_b(delta)^2
                      = sum_b 4*sinh(delta*log b)^2 / b.

   A common nonzero reflected-Green factor E multiplies every coordinate and
   therefore multiplies every energy below by E^2 without changing the locus.

3. Canonical provenance-state energy
   The centered carry defect axis of a base b has norm^2 (b-1)/b and its
   canonical dual has norm^2 b/(b-1).  Thus the exact finite state ledger is

       state_energy = sum_b [b/(b-1)] * a_b(delta)^2,

   so readout_energy <= state_energy <= 2*readout_energy for every atlas.

The script does NOT claim that every scalar zero automatically produces this
state.  It measures the finite radial/readout obstruction and its canonical
finite provenance realization.

Hybrid architecture
-------------------
* CPU: ordered, segmented, multithread prime sieve (for the prime subatlas).
* GPU: streamed CuPy float64 accumulation for all selected camera families.
* Memory: bounded segments; no global base, log, or cumulative arrays.
* Reliability: exact axis zero, CUDA OOM backoff by segment-size selection,
  checkpoint ledgers, atomic CSV/JSON writes, geometry/self tests.

Even-camera geometry
--------------------
For a natural base b >= 3 the full saturated radii are r=1,...,floor(b/2).
For even b, r=b/2 is antipodal: +r == -r modulo b, but c-r and c+r are
geometrically distinct legs.  In C4, r=1 reads odd neighbors and r=2 reads
numbers divisible by 2 but not by 4.  C2 remains the special aligned chart
with centers 4m and only the r=1 sector.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import platform
import time
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence

import numpy as np

try:
    import cupy as cp

    HAS_CUPY = True
except Exception:  # pragma: no cover - environment dependent
    cp = None
    HAS_CUPY = False


FAMILY_ORDER = ("all", "odd", "even", "prime")
DEFAULT_THRESHOLDS = ((13.0 / 250.0) ** 2, 1.0)
FLOAT_TINY = np.finfo(np.float64).tiny


@dataclass(frozen=True)
class BaseSegment:
    index: int
    low: int
    high: int
    primes: np.ndarray
    sieve_seconds: float


@dataclass
class RuntimeStats:
    segments: int = 0
    sieve_cpu_seconds: float = 0.0
    compute_seconds: float = 0.0
    processed_bases: int = 0
    processed_primes: int = 0


@dataclass
class FamilyBatch:
    count: int
    h_total: float
    l2_total: float
    readout_totals: dict[float, float]
    state_totals: dict[float, float]
    h_at: dict[int, float]
    l2_at: dict[int, float]
    readout_at: dict[float, dict[int, float]]
    state_at: dict[float, dict[int, float]]
    witness_local: dict[tuple[str, float, float], int]


def parse_float_list(text: str, name: str) -> list[float]:
    try:
        values = [float(piece.strip()) for piece in text.split(",") if piece.strip()]
    except ValueError as exc:
        raise SystemExit(f"{name} invalido: {exc}") from exc
    if not values or not all(math.isfinite(value) for value in values):
        raise SystemExit(f"{name} precisa conter numeros finitos")
    return values


def parse_int_list(text: str, name: str) -> list[int]:
    try:
        values = [int(piece.strip()) for piece in text.split(",") if piece.strip()]
    except ValueError as exc:
        raise SystemExit(f"{name} invalido: {exc}") from exc
    if not values:
        raise SystemExit(f"{name} nao pode ser vazio")
    return values


def parse_families(text: str) -> tuple[str, ...]:
    requested = [piece.strip().lower() for piece in text.split(",") if piece.strip()]
    invalid = sorted(set(requested) - set(FAMILY_ORDER))
    if invalid:
        raise SystemExit(f"familias invalidas: {invalid}; use {','.join(FAMILY_ORDER)}")
    if not requested:
        raise SystemExit("selecione pelo menos uma familia")
    return tuple(family for family in FAMILY_ORDER if family in set(requested))


def unique_abs_deltas(deltas: Sequence[float]) -> tuple[list[float], dict[float, float]]:
    unique: list[float] = []
    mapping: dict[float, float] = {}
    seen: set[float] = set()
    for raw in deltas:
        value = float(raw)
        canonical = abs(value)
        mapping[value] = canonical
        if canonical > 0.0 and canonical not in seen:
            seen.add(canonical)
            unique.append(canonical)
    return unique, mapping


def build_checkpoints(bmax: int) -> np.ndarray:
    points: set[int] = {bmax}
    for value in (10, 20, 30, 50, 70, 100, 200, 300, 500, 700):
        if 2 <= value <= bmax:
            points.add(value)
    power = 3
    while 10**power <= bmax * 10:
        scale = 10**power
        for multiplier in (1, 2, 3, 5, 7):
            value = multiplier * scale
            if 2 <= value <= bmax:
                points.add(value)
        power += 1
    return np.asarray(sorted(points), dtype=np.int64)


def simple_prime_sieve(limit: int) -> np.ndarray:
    if limit < 2:
        return np.empty(0, dtype=np.int64)
    if limit == 2:
        return np.array([2], dtype=np.int64)
    size = (limit - 1) // 2
    mask = np.ones(size, dtype=np.bool_)
    root = math.isqrt(limit)
    for value in range(3, root + 1, 2):
        index = (value - 3) // 2
        if mask[index]:
            start = (value * value - 3) // 2
            mask[start::value] = False
    odds = 2 * np.flatnonzero(mask).astype(np.int64) + 3
    return np.concatenate((np.array([2], dtype=np.int64), odds))


def sieve_segment(
    index: int,
    low: int,
    high: int,
    base_odd_primes: np.ndarray,
) -> BaseSegment:
    started = time.perf_counter()
    pieces: list[np.ndarray] = []
    if low <= 2 <= high:
        pieces.append(np.array([2], dtype=np.int64))

    odd_low = max(3, low)
    if odd_low % 2 == 0:
        odd_low += 1
    odd_high = high if high % 2 else high - 1
    if odd_low <= odd_high:
        count = (odd_high - odd_low) // 2 + 1
        mask = np.ones(count, dtype=np.bool_)
        root = math.isqrt(odd_high)
        for raw in base_odd_primes:
            prime = int(raw)
            if prime > root:
                break
            start = max(prime * prime, ((odd_low + prime - 1) // prime) * prime)
            if start % 2 == 0:
                start += prime
            mask[(start - odd_low) // 2 :: prime] = False
        pieces.append((odd_low + 2 * np.flatnonzero(mask)).astype(np.int64))

    primes = np.concatenate(pieces) if pieces else np.empty(0, dtype=np.int64)
    return BaseSegment(index, low, high, primes, time.perf_counter() - started)


def segment_bounds(bmax: int, segment_bases: int) -> list[tuple[int, int, int]]:
    output: list[tuple[int, int, int]] = []
    low = 2
    index = 0
    while low <= bmax:
        high = min(bmax, low + segment_bases - 1)
        output.append((index, low, high))
        low = high + 1
        index += 1
    return output


def iter_segments(
    bmax: int,
    segment_bases: int,
    workers: int,
    prefetch: int,
    need_primes: bool,
) -> Iterator[BaseSegment]:
    bounds = segment_bounds(bmax, segment_bases)
    if not need_primes:
        for index, low, high in bounds:
            yield BaseSegment(index, low, high, np.empty(0, dtype=np.int64), 0.0)
        return

    base = simple_prime_sieve(math.isqrt(bmax))
    base_odd = base[base != 2]
    if workers <= 1:
        for index, low, high in bounds:
            yield sieve_segment(index, low, high, base_odd)
        return

    window = max(workers, prefetch)
    futures: dict[int, Future[BaseSegment]] = {}
    next_submit = 0
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="camera-prime-sieve") as pool:
        while next_submit < min(window, len(bounds)):
            index, low, high = bounds[next_submit]
            futures[next_submit] = pool.submit(sieve_segment, index, low, high, base_odd)
            next_submit += 1
        for position in range(len(bounds)):
            result = futures.pop(position).result()
            if next_submit < len(bounds):
                index, low, high = bounds[next_submit]
                futures[next_submit] = pool.submit(sieve_segment, index, low, high, base_odd)
                next_submit += 1
            yield result


def family_values_numpy(family: str, low: int, high: int, primes: np.ndarray) -> np.ndarray:
    if family == "prime":
        return np.asarray(primes, dtype=np.int64)
    values = np.arange(low, high + 1, dtype=np.int64)
    if family == "odd":
        return values[(values & 1) == 1]
    if family == "even":
        return values[(values & 1) == 0]
    return values


def selected_cumulative(cumulative, indices: np.ndarray, xp) -> dict[int, float]:
    if indices.size == 0:
        return {}
    if xp is np:
        values = np.asarray(cumulative[indices], dtype=np.float64)
    else:
        values = cp.asnumpy(cumulative[cp.asarray(indices, dtype=cp.int64)])
    return {int(index): float(value) for index, value in zip(indices, values, strict=True)}


def process_family(
    family: str,
    values_np: np.ndarray,
    abs_deltas: Sequence[float],
    take_indices: np.ndarray,
    readout_before: dict[float, float],
    state_before: dict[float, float],
    unresolved: Sequence[tuple[str, float, float]],
    green_energy_sq: float,
    xp,
) -> FamilyBatch:
    if values_np.size == 0:
        return FamilyBatch(
            count=0,
            h_total=0.0,
            l2_total=0.0,
            readout_totals={delta: 0.0 for delta in abs_deltas},
            state_totals={delta: 0.0 for delta in abs_deltas},
            h_at={},
            l2_at={},
            readout_at={delta: {} for delta in abs_deltas},
            state_at={delta: {} for delta in abs_deltas},
            witness_local={},
        )

    values = xp.asarray(values_np, dtype=xp.float64)
    logs = xp.log(values)
    inv = 1.0 / values
    h_cum = xp.cumsum(inv, dtype=xp.float64)
    l2_cum = xp.cumsum(logs * logs * inv, dtype=xp.float64)

    h_total = float(h_cum[-1]) if xp is np else float(h_cum[-1].get())
    l2_total = float(l2_cum[-1]) if xp is np else float(l2_cum[-1].get())
    h_at = selected_cumulative(h_cum, take_indices, xp)
    l2_at = selected_cumulative(l2_cum, take_indices, xp)

    unresolved_set = {(f, d, t) for f, d, t in unresolved if f == family}
    readout_totals: dict[float, float] = {}
    state_totals: dict[float, float] = {}
    readout_at: dict[float, dict[int, float]] = {}
    state_at: dict[float, dict[int, float]] = {}
    witness_local: dict[tuple[str, float, float], int] = {}

    dual_norm_sq = values / (values - 1.0)
    for delta in abs_deltas:
        x = np.float64(delta) * logs
        radial = 2.0 * xp.sinh(x)
        term = green_energy_sq * radial * radial * inv
        readout_cum = xp.cumsum(term, dtype=xp.float64)
        state_cum = xp.cumsum(term * dual_norm_sq, dtype=xp.float64)

        readout_total = (
            float(readout_cum[-1]) if xp is np else float(readout_cum[-1].get())
        )
        state_total = float(state_cum[-1]) if xp is np else float(state_cum[-1].get())
        readout_totals[delta] = readout_total
        state_totals[delta] = state_total
        readout_at[delta] = selected_cumulative(readout_cum, take_indices, xp)
        state_at[delta] = selected_cumulative(state_cum, take_indices, xp)

        before = readout_before[delta]
        after = before + readout_total
        for key in unresolved_set:
            _, key_delta, threshold = key
            if key_delta != delta or not (before < threshold <= after):
                continue
            target = np.float64(threshold - before)
            if xp is np:
                local = int(np.searchsorted(readout_cum, target, side="left"))
            else:
                local = int(cp.searchsorted(readout_cum, target, side="left").get())
            witness_local[key] = local

        del x, radial, term, readout_cum, state_cum

    if xp is not np:
        cp.cuda.Stream.null.synchronize()
    return FamilyBatch(
        count=int(values_np.size),
        h_total=h_total,
        l2_total=l2_total,
        readout_totals=readout_totals,
        state_totals=state_totals,
        h_at=h_at,
        l2_at=l2_at,
        readout_at=readout_at,
        state_at=state_at,
        witness_local=witness_local,
    )


class AtlasLedger:
    def __init__(
        self,
        families: Sequence[str],
        deltas: Sequence[float],
        checkpoints: np.ndarray,
    ) -> None:
        self.families = tuple(families)
        self.deltas = [float(value) for value in deltas]
        self.abs_deltas, self.delta_map = unique_abs_deltas(self.deltas)
        self.checkpoints = checkpoints
        n = checkpoints.size

        self.count = {family: 0 for family in self.families}
        self.h = {family: 0.0 for family in self.families}
        self.l2 = {family: 0.0 for family in self.families}
        self.readout = {
            family: {delta: 0.0 for delta in self.abs_deltas}
            for family in self.families
        }
        self.state = {
            family: {delta: 0.0 for delta in self.abs_deltas}
            for family in self.families
        }

        self.count_at = {family: np.zeros(n, dtype=np.int64) for family in self.families}
        self.h_at = {family: np.zeros(n, dtype=np.float64) for family in self.families}
        self.l2_at = {family: np.zeros(n, dtype=np.float64) for family in self.families}
        self.readout_at = {
            family: {delta: np.zeros(n, dtype=np.float64) for delta in self.abs_deltas}
            for family in self.families
        }
        self.state_at = {
            family: {delta: np.zeros(n, dtype=np.float64) for delta in self.abs_deltas}
            for family in self.families
        }
        self.done = np.zeros(n, dtype=np.bool_)
        self.witnesses: dict[tuple[str, float, float], dict[str, int] | None] = {
            (family, delta, threshold): None
            for family in self.families
            for delta in self.abs_deltas
            for threshold in DEFAULT_THRESHOLDS
        }

    def unresolved(self) -> list[tuple[str, float, float]]:
        return [key for key, value in self.witnesses.items() if value is None]

    def consume_segment(
        self,
        segment: BaseSegment,
        green_energy_sq: float,
        xp,
    ) -> None:
        checkpoint_positions = np.flatnonzero(
            (~self.done)
            & (self.checkpoints >= segment.low)
            & (self.checkpoints <= segment.high)
        )
        unresolved = self.unresolved()

        family_results: dict[str, tuple[np.ndarray, FamilyBatch, dict[int, int]]] = {}
        for family in self.families:
            values = family_values_numpy(family, segment.low, segment.high, segment.primes)
            local_by_position: dict[int, int] = {}
            requested: list[int] = []
            for position in checkpoint_positions:
                checkpoint = int(self.checkpoints[position])
                local = int(np.searchsorted(values, checkpoint, side="right") - 1)
                local_by_position[int(position)] = local
                if local >= 0:
                    requested.append(local)
            take_indices = np.unique(np.asarray(requested, dtype=np.int64))
            result = process_family(
                family,
                values,
                self.abs_deltas,
                take_indices,
                self.readout[family],
                self.state[family],
                unresolved,
                green_energy_sq,
                xp,
            )
            family_results[family] = (values, result, local_by_position)

        for position in checkpoint_positions:
            for family in self.families:
                values, result, local_by_position = family_results[family]
                local = local_by_position[int(position)]
                if local < 0:
                    self.count_at[family][position] = self.count[family]
                    self.h_at[family][position] = self.h[family]
                    self.l2_at[family][position] = self.l2[family]
                    for delta in self.abs_deltas:
                        self.readout_at[family][delta][position] = self.readout[family][delta]
                        self.state_at[family][delta][position] = self.state[family][delta]
                else:
                    self.count_at[family][position] = self.count[family] + local + 1
                    self.h_at[family][position] = self.h[family] + result.h_at[local]
                    self.l2_at[family][position] = self.l2[family] + result.l2_at[local]
                    for delta in self.abs_deltas:
                        self.readout_at[family][delta][position] = (
                            self.readout[family][delta] + result.readout_at[delta][local]
                        )
                        self.state_at[family][delta][position] = (
                            self.state[family][delta] + result.state_at[delta][local]
                        )
            self.done[position] = True

        for family in self.families:
            values, result, _ = family_results[family]
            count_before = self.count[family]
            for key, local in result.witness_local.items():
                if self.witnesses[key] is None:
                    self.witnesses[key] = {
                        "member_count": count_before + local + 1,
                        "base": int(values[local]),
                    }
            self.count[family] += result.count
            self.h[family] += result.h_total
            self.l2[family] += result.l2_total
            for delta in self.abs_deltas:
                self.readout[family][delta] += result.readout_totals[delta]
                self.state[family][delta] += result.state_totals[delta]

    def finalize(self) -> None:
        for position in np.flatnonzero(~self.done):
            for family in self.families:
                self.count_at[family][position] = self.count[family]
                self.h_at[family][position] = self.h[family]
                self.l2_at[family][position] = self.l2[family]
                for delta in self.abs_deltas:
                    self.readout_at[family][delta][position] = self.readout[family][delta]
                    self.state_at[family][delta][position] = self.state[family][delta]
            self.done[position] = True


def cuda_info(device: int) -> dict[str, object]:
    assert HAS_CUPY and cp is not None
    count = int(cp.cuda.runtime.getDeviceCount())
    if not 0 <= device < count:
        raise RuntimeError(f"CUDA device {device} invalido; detectados {count}")
    cp.cuda.Device(device).use()
    properties = cp.cuda.runtime.getDeviceProperties(device)
    raw_name = properties.get("name", b"CUDA GPU")
    name = raw_name.decode(errors="replace") if isinstance(raw_name, bytes) else str(raw_name)
    free_bytes, total_bytes = cp.cuda.runtime.memGetInfo()
    return {
        "device": device,
        "name": name,
        "free_bytes": int(free_bytes),
        "total_bytes": int(total_bytes),
        "cupy": cp.__version__,
    }


def resolve_backend(requested: str, device: int) -> tuple[object, str, dict[str, object]]:
    if requested == "cpu":
        return np, "CPU", {}
    if HAS_CUPY and cp is not None:
        try:
            info = cuda_info(device)
            return cp, "CUDA", info
        except Exception as exc:
            if requested == "cuda":
                raise SystemExit(f"CUDA indisponivel: {exc}") from exc
    elif requested == "cuda":
        raise SystemExit("CUDA pedido, mas CuPy nao esta instalado")
    return np, "CPU(auto)", {}


def auto_segment_bases(xp, cuda_meta: dict[str, object], requested: int) -> int:
    if requested > 0:
        return requested
    if xp is np:
        return 1_000_000
    free_bytes = int(cuda_meta["free_bytes"])
    # Per family and per delta arrays are released before the next family/delta.
    # 96 bytes/base is deliberately conservative for values/log/inv/cumsums.
    estimate = int((0.35 * free_bytes) // 96)
    return max(250_000, min(5_000_000, estimate))


def camera_geometry(base: int) -> dict[str, object]:
    if base < 2:
        raise ValueError("base precisa ser >= 2")
    if base == 2:
        return {
            "base": 2,
            "kind": "C2_special_aligned",
            "center_step": 4,
            "radii": [1],
            "antipodal_radius": None,
            "description": "centros 4m; somente setor impar r=1",
        }
    half = base // 2
    antipodal = half if base % 2 == 0 else None
    return {
        "base": base,
        "kind": "natural_saturated_even" if antipodal else "natural_saturated_odd",
        "center_step": base,
        "radii": list(range(1, half + 1)),
        "antipodal_radius": antipodal,
        "description": (
            f"centros {base}m; raios 1..{half}; r={half} antipodal"
            if antipodal
            else f"centros {base}m; raios 1..{half}"
        ),
    }


def print_geometry(bases: Sequence[int]) -> None:
    print(" Camera geometry audit:")
    for base in bases:
        row = camera_geometry(base)
        suffix = ""
        if base == 4:
            suffix = "; r=1 -> impares, r=2 -> pares com v2=1"
        print(
            f"   C{base:<3d} {row['kind']:<28s} "
            f"centers={row['center_step']}m radii={row['radii']}{suffix}"
        )


def direct_reference(
    bmax: int,
    family: str,
    delta: float,
    green_energy_sq: float,
) -> tuple[float, float]:
    bases = np.arange(2, bmax + 1, dtype=np.int64)
    if family == "odd":
        bases = bases[bases % 2 == 1]
    elif family == "even":
        bases = bases[bases % 2 == 0]
    elif family == "prime":
        primes = simple_prime_sieve(bmax)
        bases = primes
    values = bases.astype(np.float64)
    radial = 2.0 * np.sinh(abs(delta) * np.log(values))
    term = green_energy_sq * radial * radial / values
    return float(np.sum(term)), float(np.sum(term * values / (values - 1.0)))


def self_test() -> int:
    c4 = camera_geometry(4)
    assert c4["radii"] == [1, 2]
    centers = 4 * np.arange(1, 20, dtype=np.int64)
    r2 = np.concatenate((centers - 2, centers + 2))
    assert np.all(r2 % 2 == 0) and np.all(r2 % 4 != 0)

    bmax = 10_000
    deltas = [0.0, 0.01, 0.2]
    families = FAMILY_ORDER
    checkpoints = np.array([100, 1000, bmax], dtype=np.int64)
    ledger = AtlasLedger(families, deltas, checkpoints)
    for segment in iter_segments(
        bmax, segment_bases=777, workers=2, prefetch=4, need_primes=True
    ):
        ledger.consume_segment(segment, green_energy_sq=1.0, xp=np)
    ledger.finalize()

    for family in families:
        assert ledger.readout_at[family].get(0.0) is None  # zero is stored implicitly
        for delta in (0.01, 0.2):
            readout, state = direct_reference(bmax, family, delta, 1.0)
            got_r = ledger.readout[family][delta]
            got_s = ledger.state[family][delta]
            if not math.isclose(got_r, readout, rel_tol=3e-13, abs_tol=1e-14):
                raise AssertionError((family, delta, got_r, readout))
            if not math.isclose(got_s, state, rel_tol=3e-13, abs_tol=1e-14):
                raise AssertionError((family, delta, got_s, state))
            assert got_r <= got_s <= 2.0 * got_r + 1e-12

    print("SELF-TEST ALL-BASE CARRY BESSEL ATLAS")
    print("  PASS  C4 has r=1 and antipodal r=2")
    print("  PASS  r=2 legs are even and not C4 centers")
    print("  PASS  segmented multithread prime sieve")
    print("  PASS  streamed family ledgers match direct NumPy")
    print("  PASS  readout <= provenance state <= 2*readout")
    return 0


def atomic_csv(path: Path, ledger: AtlasLedger) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    fields = ["cutoff"]
    for family in ledger.families:
        fields.extend((f"{family}:count", f"{family}:H", f"{family}:L2"))
        for raw in ledger.deltas:
            fields.extend((f"{family}:readout:{raw}", f"{family}:state:{raw}"))
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for position, cutoff in enumerate(ledger.checkpoints):
            row: dict[str, int | float] = {"cutoff": int(cutoff)}
            for family in ledger.families:
                row[f"{family}:count"] = int(ledger.count_at[family][position])
                row[f"{family}:H"] = float(ledger.h_at[family][position])
                row[f"{family}:L2"] = float(ledger.l2_at[family][position])
                for raw in ledger.deltas:
                    canonical = ledger.delta_map[raw]
                    if canonical == 0.0:
                        readout = state = 0.0
                    else:
                        readout = ledger.readout_at[family][canonical][position]
                        state = ledger.state_at[family][canonical][position]
                    row[f"{family}:readout:{raw}"] = float(readout)
                    row[f"{family}:state:{raw}"] = float(state)
            writer.writerow(row)
    temporary.replace(path)


def atomic_json(
    path: Path,
    ledger: AtlasLedger,
    metadata: dict[str, object],
    geometries: Sequence[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    rows = []
    for position, cutoff in enumerate(ledger.checkpoints):
        item: dict[str, object] = {"cutoff": int(cutoff), "families": {}}
        family_payload: dict[str, object] = {}
        for family in ledger.families:
            energies: dict[str, object] = {}
            for raw in ledger.deltas:
                canonical = ledger.delta_map[raw]
                if canonical == 0.0:
                    readout = state = 0.0
                else:
                    readout = float(ledger.readout_at[family][canonical][position])
                    state = float(ledger.state_at[family][canonical][position])
                energies[str(raw)] = {"readout": readout, "state": state}
            family_payload[family] = {
                "count": int(ledger.count_at[family][position]),
                "H": float(ledger.h_at[family][position]),
                "L2": float(ledger.l2_at[family][position]),
                "energies": energies,
            }
        item["families"] = family_payload
        rows.append(item)

    witnesses = []
    for key, value in ledger.witnesses.items():
        family, delta, threshold = key
        witnesses.append(
            {
                "family": family,
                "abs_delta": delta,
                "threshold": threshold,
                "witness": value,
            }
        )
    payload = {
        "metadata": metadata,
        "camera_geometries": list(geometries),
        "rows": rows,
        "witnesses": witnesses,
    }
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="All-base carry Green/Bessel radial atlas: CPU sieve + CUDA float64",
    )
    parser.add_argument("--bmax", type=float, default=1e8, help="largest positional base")
    parser.add_argument(
        "--deltas", default="0,0.001,0.01,0.05,0.1,0.25,0.4",
        help="critical displacements sigma-1/2",
    )
    parser.add_argument(
        "--families", default="all,odd,even,prime",
        help="camera families: all,odd,even,prime",
    )
    parser.add_argument("--backend", choices=("auto", "cpu", "cuda"), default="auto")
    parser.add_argument("--cuda-device", type=int, default=0)
    parser.add_argument(
        "--sieve-workers", type=int, default=max(1, (os.cpu_count() or 1) - 1),
        help="CPU threads for the prime reference sieve",
    )
    parser.add_argument("--sieve-prefetch", type=int, default=0)
    parser.add_argument(
        "--segment-bases", type=int, default=0,
        help="bases per streamed segment; 0 chooses from RAM/VRAM",
    )
    parser.add_argument(
        "--green-energy", type=float, default=1.0,
        help="common reflected-Green amplitude E; energies are multiplied by E^2",
    )
    parser.add_argument("--geometry-bases", default="2,3,4,5,6,8,9")
    parser.add_argument("--progress-seconds", type=float, default=2.0)
    parser.add_argument("--validate-sieve", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--csv-out", type=Path)
    parser.add_argument("--json-out", type=Path)
    return parser


def validate_sieve(bmax: int, workers: int, prefetch: int) -> None:
    limit = min(bmax, 5_000_000)
    pieces = [
        segment.primes
        for segment in iter_segments(limit, 250_000, workers, prefetch, True)
        if segment.primes.size
    ]
    got = np.concatenate(pieces) if pieces else np.empty(0, dtype=np.int64)
    expected = simple_prime_sieve(limit)
    if not np.array_equal(got, expected):
        raise RuntimeError(f"crivo segmentado divergiu ate {limit}")
    print(f" sieve validation    : PASS through {limit:,}")


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.self_test:
        return self_test()

    bmax = int(args.bmax)
    if bmax < 2:
        raise SystemExit("bmax precisa ser >= 2")
    deltas = parse_float_list(args.deltas, "deltas")
    families = parse_families(args.families)
    geometry_bases = parse_int_list(args.geometry_bases, "geometry-bases")
    if any(base < 2 for base in geometry_bases):
        raise SystemExit("toda base geometrica precisa ser >= 2")
    if not math.isfinite(args.green_energy):
        raise SystemExit("green-energy precisa ser finito")
    green_energy_sq = float(args.green_energy) ** 2

    max_x = max(abs(delta) for delta in deltas) * math.log(float(bmax))
    if max_x > 350.0:
        raise SystemExit(
            f"abs(delta)*log(bmax)={max_x:.2f} excede a faixa segura float64 para sinh^2"
        )

    xp, backend_name, cuda_meta = resolve_backend(args.backend, args.cuda_device)
    segment_bases = auto_segment_bases(xp, cuda_meta, args.segment_bases)
    workers = max(1, int(args.sieve_workers))
    prefetch = args.sieve_prefetch if args.sieve_prefetch > 0 else 2 * workers
    checkpoints = build_checkpoints(bmax)
    ledger = AtlasLedger(families, deltas, checkpoints)
    need_primes = "prime" in families

    print("=" * 116)
    print(" ALL-BASE CARRY GREEN/BESSEL ATLAS — MULTITHREAD CPU SIEVE + STREAMED CUDA")
    print("=" * 116)
    print(f" backend             : {backend_name}")
    if cuda_meta:
        print(f" CUDA device         : {cuda_meta['device']} — {cuda_meta['name']}")
        print(
            f" CUDA memory         : {int(cuda_meta['free_bytes']) / 2**30:.2f} GiB free / "
            f"{int(cuda_meta['total_bytes']) / 2**30:.2f} GiB total"
        )
    print(f" base cutoff         : 2 <= b <= {bmax:,}")
    print(f" camera families     : {families}")
    print(f" deltas              : {deltas}")
    print(f" common Green E      : {args.green_energy:g}  (energy factor E^2={green_energy_sq:g})")
    print(f" segment size        : {segment_bases:,} bases")
    print(f" prime sieve threads : {workers if need_primes else 0}")
    print(f" sieve prefetch      : {prefetch if need_primes else 0}")
    print(" radial coordinate   : b^(-1/2) * (b^delta - b^(-delta))")
    print(" state ledger        : sum [b/(b-1)] * readout_coordinate^2")
    print(" axis guarantee      : readout <= state <= 2*readout")
    print("-")
    print_geometry(geometry_bases)
    print("-" * 116)

    if args.validate_sieve and need_primes:
        validate_sieve(bmax, workers, prefetch)

    stats = RuntimeStats()
    started = time.perf_counter()
    last_progress = started
    for segment in iter_segments(bmax, segment_bases, workers, prefetch, need_primes):
        stats.segments += 1
        stats.sieve_cpu_seconds += segment.sieve_seconds
        stats.processed_bases += segment.high - segment.low + 1
        stats.processed_primes += int(segment.primes.size)
        compute_started = time.perf_counter()
        ledger.consume_segment(segment, green_energy_sq, xp)
        stats.compute_seconds += time.perf_counter() - compute_started

        now = time.perf_counter()
        if args.progress_seconds > 0 and now - last_progress >= args.progress_seconds:
            elapsed = now - started
            print(
                f"[progress] {100.0 * segment.high / bmax:6.2f}% "
                f"b<={segment.high:,} segments={stats.segments:,} "
                f"rate={stats.processed_bases / elapsed / 1e6:.2f} Mbase/s",
                flush=True,
            )
            last_progress = now
    ledger.finalize()
    wall = time.perf_counter() - started

    print("-" * 116)
    print(f" wall elapsed        : {wall:.3f}s")
    print(f" compute wall        : {stats.compute_seconds:.3f}s")
    print(f" sieve CPU work      : {stats.sieve_cpu_seconds:.3f}s (sum of worker tasks)")
    print(f" streamed segments   : {stats.segments:,}")
    print(f" bases processed     : {stats.processed_bases:,}")
    if need_primes:
        print(f" primes classified   : {stats.processed_primes:,}")
    if wall > 0:
        print(f" throughput          : {stats.processed_bases / wall / 1e6:.3f} million bases/s")

    for family in families:
        print("-" * 116)
        print(f" FAMILY {family.upper()} — final cutoff {bmax:,}, members={ledger.count[family]:,}")
        print(f" {'delta':>10} {'readout energy':>20} {'state energy':>20} {'state/readout':>16} {'linearized ratio':>18}")
        for raw in deltas:
            canonical = ledger.delta_map[raw]
            if canonical == 0.0:
                readout = state = 0.0
                state_ratio = math.nan
                linear_ratio = math.nan
            else:
                readout = ledger.readout[family][canonical]
                state = ledger.state[family][canonical]
                state_ratio = state / readout if readout > 0.0 else math.nan
                linear_prediction = 4.0 * canonical * canonical * green_energy_sq * ledger.l2[family]
                linear_ratio = readout / linear_prediction if linear_prediction > 0.0 else math.nan
            print(
                f" {raw:>10g} {readout:>20.10e} {state:>20.10e} "
                f"{state_ratio:>16.9f} {linear_ratio:>18.9f}"
            )
        print(
            "   linearized ratio = readout / (4 delta^2 E^2 sum(log(b)^2/b)); "
            "near 1 only when |delta| log(bmax) is small."
        )

    print("-" * 116)
    print(" First radial-readout threshold witnesses:")
    for family in families:
        for raw in deltas:
            canonical = ledger.delta_map[raw]
            if canonical == 0.0:
                continue
            pieces = []
            for threshold in DEFAULT_THRESHOLDS:
                witness = ledger.witnesses[(family, canonical, threshold)]
                pieces.append(
                    f"{threshold:.6g}->"
                    + (
                        f"base {witness['base']:,} (member #{witness['member_count']:,})"
                        if witness is not None
                        else ">bmax"
                    )
                )
            print(f"   {family:>5s} delta={raw:g}: " + ", ".join(pieces))

    geometries = [camera_geometry(base) for base in geometry_bases]
    metadata: dict[str, object] = {
        "schema": "org.native-carry.all-base-green-bessel-radial/v2",
        "bmax": bmax,
        "families": list(families),
        "deltas": deltas,
        "green_energy": args.green_energy,
        "backend": backend_name,
        "segment_bases": segment_bases,
        "sieve_workers": workers if need_primes else 0,
        "wall_seconds": wall,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "cuda": cuda_meta or None,
        "prime_theorem_scope": "prime family",
        "all_base_status": "finite numerical/generalized ledger; not asserted as an existing Lean theorem",
    }
    if args.csv_out:
        atomic_csv(args.csv_out, ledger)
        print(f" CSV written         : {args.csv_out}")
    if args.json_out:
        atomic_json(args.json_out, ledger, metadata, geometries)
        print(f" JSON written        : {args.json_out}")

    print("=" * 116)
    print(" VERDICT")
    print(" * C4 is a full even camera with r=1 and antipodal r=2; it is not C2.")
    print(" * The prime series is a subatlas, not the definition of the positional atlas.")
    print(" * The program measures radial Green readout and finite provenance-state energy, not the missing zero-to-state handoff.")
    print("=" * 116)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
