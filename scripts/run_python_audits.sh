#!/usr/bin/env bash
set -euo pipefail

export PYTHONDONTWRITEBYTECODE=1

python3 laboratory/native_carry_primitive_real_operator_all_bases.py --self-test
python3 laboratory/native_carry_precision_ladder.py --self-test
python3 laboratory/native_carry_asymptotic_corrector.py --self-test
python3 scripts/verify_interval_certificate.py \
  certificates/c3/c3_m16364_raw_minimum.json --recompute
python3 scripts/verify_c3_tail_certificate.py \
  certificates/c3/c3_m16364_tail_limit_minimum.json --recompute
python3 scripts/verify_c3_oriented_tail_certificate.py \
  certificates/c3/c3_m16364_oriented_limit_minimum.json --recompute
python3 scripts/verify_c3_oriented_tail_certificate.py \
  certificates/c3/c3_m65536_oriented_limit_minimum.json --recompute
python3 scripts/verify_c3_oriented_tail_certificate.py \
  certificates/c3/c3_m131072_oriented_limit_minimum.json --recompute
python3 scripts/verify_c3_contraction_ladder.py \
  certificates/c3/c3_oriented_contraction_ladder.json --recompute
python3 scripts/verify_c3_uniform_residual.py \
  certificates/c3/c3_uniform_residual_decomposition.json --recompute
python3 scripts/verify_c3_stationary_localization.py \
  certificates/c3/c3_stationary_localization.json --recompute
sha256sum --check audit/SOURCE_SHA256.txt
sha256sum --check audit/CERTIFICATE_SHA256.txt
python3 -m pytest -q
