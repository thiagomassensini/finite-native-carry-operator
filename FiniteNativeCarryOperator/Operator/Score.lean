import FiniteNativeCarryOperator.Operator.FiniteReal
import Mathlib.Tactic

namespace FiniteNativeCarryOperator.Operator

noncomputable section

/-- Euclidean quadratic energy in the real carrier plane. -/
def quadraticEnergy (u : RealPlane) : ℝ := u.1 ^ 2 + u.2 ^ 2

/-- Euclidean scalar product in the native real plane. -/
def realDot (u v : RealPlane) : ℝ := u.1 * v.1 + u.2 * v.2

/-- Oriented real-plane determinant. -/
def realDet (u v : RealPlane) : ℝ := u.1 * v.2 - u.2 * v.1

/-- The published dimensionless score, parameterized by its positive ledger terms. -/
def normalizedScore (coordinateCount totalEnergy : ℝ) (u : RealPlane) : ℝ :=
  quadraticEnergy u / (coordinateCount * totalEnergy)

/-- Canonical finite zero predicate for the already weighted operator. -/
def IsFiniteNativeCarryOperatorZero
    (camera cutoff : ℕ) (time : ℝ) : Prop :=
  finiteNativeOperator camera cutoff time = 0

theorem quadraticEnergy_nonneg (u : RealPlane) : 0 ≤ quadraticEnergy u := by
  exact add_nonneg (sq_nonneg u.1) (sq_nonneg u.2)

theorem quadraticEnergy_eq_zero_iff (u : RealPlane) :
    quadraticEnergy u = 0 ↔ u = 0 := by
  rcases u with ⟨x, y⟩
  change x ^ 2 + y ^ 2 = 0 ↔ (x, y) = (0, 0)
  constructor
  · intro h
    have hx : x = 0 := by nlinarith [sq_nonneg x, sq_nonneg y]
    have hy : y = 0 := by nlinarith [sq_nonneg x, sq_nonneg y]
    simp [hx, hy]
  · intro h
    have hx : x = 0 := congrArg Prod.fst h
    have hy : y = 0 := congrArg Prod.snd h
    simp [hx, hy]

/-- Real two-dimensional Lagrange identity. -/
theorem realDot_sq_add_realDet_sq (u v : RealPlane) :
    realDot u v ^ 2 + realDet u v ^ 2 =
      quadraticEnergy u * quadraticEnergy v := by
  rcases u with ⟨ux, uy⟩
  rcases v with ⟨vx, vy⟩
  simp only [realDot, realDet, quadraticEnergy]
  ring

/-- Dot and determinant against one nonzero vector determine the other vector. -/
theorem eq_zero_iff_realDot_eq_zero_and_realDet_eq_zero_of_right_ne_zero
    (u v : RealPlane) (hv : v ≠ 0) :
    u = 0 ↔ realDot u v = 0 ∧ realDet u v = 0 := by
  constructor
  · intro hu
    simp [hu, realDot, realDet]
  · rintro ⟨hdot, hdet⟩
    have hvEnergy : quadraticEnergy v ≠ 0 := by
      intro hzero
      exact hv ((quadraticEnergy_eq_zero_iff v).mp hzero)
    have hproduct : quadraticEnergy u * quadraticEnergy v = 0 := by
      rw [← realDot_sq_add_realDet_sq u v, hdot, hdet]
      norm_num
    have huEnergy : quadraticEnergy u = 0 :=
      (mul_eq_zero.mp hproduct).resolve_right hvEnergy
    exact (quadraticEnergy_eq_zero_iff u).mp huEnergy

/-- At stationarity, one determinant decides the vector zero when velocity is nonzero. -/
theorem eq_zero_iff_realDet_eq_zero_of_realDot_eq_zero_of_right_ne_zero
    (u v : RealPlane) (hv : v ≠ 0) (hdot : realDot u v = 0) :
    u = 0 ↔ realDet u v = 0 := by
  rw [eq_zero_iff_realDot_eq_zero_and_realDet_eq_zero_of_right_ne_zero u v hv]
  simp [hdot]

theorem normalizedScore_eq_zero_iff
    (coordinateCount totalEnergy : ℝ)
    (hcount : coordinateCount ≠ 0) (henergy : totalEnergy ≠ 0)
    (u : RealPlane) :
    normalizedScore coordinateCount totalEnergy u = 0 ↔ u = 0 := by
  rw [normalizedScore, div_eq_zero_iff]
  simp [hcount, henergy, quadraticEnergy_eq_zero_iff]

theorem finite_energy_eq_zero_iff
    (camera cutoff : ℕ) (time : ℝ) :
    quadraticEnergy (finiteNativeOperator camera cutoff time) = 0 ↔
      IsFiniteNativeCarryOperatorZero camera cutoff time := by
  simpa [IsFiniteNativeCarryOperatorZero] using
    quadraticEnergy_eq_zero_iff (finiteNativeOperator camera cutoff time)

end
end FiniteNativeCarryOperator.Operator
